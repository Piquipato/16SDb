import typing as tp
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import numpy as np
import pandas as pd

def extract_metadata(fragment: SeqRecord) -> tp.Dict[str, str]:
    props = {"id": "","name": "","description": ""}
    for key in props.keys():
        props[key] = getattr(fragment, key)
    props["region"] = re.search(r"Region: ([\w-]+)", props["description"]).group(1)
    props["ggid"] = re.search(r"ID: (\w+)", props["description"]).group(1)
    props["ggname"] = re.search(r"Name: (\w+)", props["description"]).group(1)
    props["ggdesc"] = re.search(
        r"Description: (.+)",
        props["description"]
    ).group(1)
    props["otu"] = re.search(r"otu_(\w+)", props["ggdesc"]).group(0)
    props["tax"] = " ".join([
        m for m in re.findall(
            r" (\w{1}__[\w ]+;|Unclassified;)",
            props["ggdesc"]
        ) if m is not None
    ])
    props["ggtag"] = props["ggdesc"].split(" ")[1]
    return props


class FragmentGenerator:

    def __init__(
            self,
            fragment_library: SeqIO,
            ratios: tp.Iterable = [],
            debug: bool = False
    ):
        # Store information on fragments
        self.fraglib = fragment_library
        self.regions = []
        self.otus_idx = []
        self.taxamap = {}
        self.metadata = []
        for f in fragment_library:
            props = extract_metadata(f)
            self.metadata.append(props)
            self.regions.append(props["region"])
            self.otus_idx.append(props["otu"])
            if not props["otu"] in self.taxamap.keys():
                self.taxamap[props["otu"]] = [props["tax"]]
            else:
                self.taxamap[props["otu"]].append(props["tax"])
        idxs = zip(self.otus_idx, self.regions)
        self.regions = np.unique(self.regions)
        self.otus = np.unique(self.otus_idx)
        self.present_matrix = pd.DataFrame(np.zeros(
            shape=(len(self.otus), len(self.regions)),
            dtype=np.int64
        ), index=self.otus, columns=self.regions)
        for i in idxs:
            self.present_matrix.loc[i]=1
        for key in self.taxamap.keys():
            self.taxamap[key] = np.unique(self.taxamap[key])
        if len(ratios) > 0:
            self.set_ratios(ratios)
        if debug:
            print(f"Regions: {self.regions}")
            print(f"Nu. OTUs: {len(self.otus)}")
            taxamap_degen = [
                key for key, val in self.taxamap.items()
                if len(val) != 1
            ]
            for ot in taxamap_degen:
                print(f"OTU: {ot}")
                print(f"Taxas: {len(self.taxamap[ot])}")
                for tax in self.taxamap[ot]:
                    print(f"\t- {tax}")
            print("Record Matrix")
            print(self.present_matrix)
            print(np.sum(self.present_matrix, axis=None) / np.prod(self.present_matrix.shape))
            if len(ratios) > 0:
                print("Ratios Matrix")
                print(self.ratios_mtx)
        
    def set_ratios(self, ratios):
        shp = ratios.shape
        if not len(shp) in [1, 2]:
            raise ValueError(f"Ratios of wrong dimensions {len(shp)}.")
        if len(shp) == 1:
            if shp[0] != len(self.regions):
                raise ValueError(
                    f"Inputed ratios {shp[0]}" \
                    f"do not match with setup ratios {len(self.ratios)}."
                )
            self.ratio_ws = pd.DataFrame(
                [ratios for _ in range(len(self.otus))],
                index=self.otus,
                columns=self.regions
            )
        if len(shp) == 2:
            if shp != (len(self.otus), len(self.regions)):
                raise ValueError(
                    f"Inputed ratios {shp}" \
                    "do not match with setup" \
                    f"ratios {self.present_matrix.shape}."
                )
            self.ratio_ws = pd.DataFrame(
                ratios,
                index=self.otus,
                columns=self.regions
            )
        self.ratio_mtx = np.multiply(self.present_matrix, self.ratio_ws)
        for otu in self.ratio_mtx:
            otu /= np.sum(otu)
        




    def generate():
        pass

    def store():
        pass



if __name__ == "__main__":
    import json
    
    fragment_library = SeqIO.parse(open("fragments.fasta"), 'fasta')
    generator = FragmentGenerator(fragment_library, debug=True)