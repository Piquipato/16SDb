import typing as tp
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scipy.special import softmax
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

def extract_metadata_library(library):
    metadata = []
    len_library = len(library)
    i = 0
    for f in SeqIO.parse(library, "fasta"):
        props = extract_metadata(f)
        props["idx"] = i
        metadata.append(props)
        i += 1
    return metadata, len_library

class FragmentGenerator:

    def __init__(
            self,
            fragment_library: str,
            ratios: tp.Iterable = [],
            debug: bool = False
    ):
        self.debug = debug
        if debug:
            print("Instantiating FragmentGenerator...")
        # Store information on fragments
        self.ratios_are_set = False
        self.fraglib = fragment_library
        self.regions = []
        self.otus_idx = []
        self.taxamap = {}
        self.metadata, self.fraglib_len = extract_metadata_library(
            fragment_library
        )
        for props in self.metadata:
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
            print("Finished loading metadata")
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
            print(np.sum(self.present_matrix.to_numpy().flatten(), axis=None) / np.prod(self.present_matrix.shape))
            if len(ratios) > 0:
                print("Ratios Matrix")
                print(self.ratios_mtx)
    
    def get_present_matrix(self):
        return self.present_matrix
    
    def get_ratios(self):
        return self.ratio_ws
    
    def get_otus(self):
        return self.otus
    
    def get_regions(self):
        return self.regions
        
    def set_ratios(self, ratios):
        if self.debug:
            print("Setting ratios...")
        shp = ratios.shape
        if not len(shp) in [1, 2]:
            raise ValueError(f"Ratios of wrong dimensions {len(shp)}.")
        if len(shp) == 1:
            if shp[0] != len(self.otus):
                raise ValueError(
                    f"Inputed ratios {shp[0]}" \
                    f"do not match with setup ratios {len(self.otus)}."
                )
            self.ratio_ws = pd.DataFrame(
                np.transpose(
                    np.array([ratios for _ in range(len(self.regions))])
                ),
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
                np.transpose(ratios),
                index=self.otus,
                columns=self.regions
            )
        self.ratio_mtx = np.multiply(self.present_matrix, self.ratio_ws)
        for col in self.ratio_mtx.columns:
            self.ratio_mtx[col] /= self.ratio_mtx[col].sum()
        self.otu_prob = pd.DataFrame(
            softmax(np.sum(self.ratio_mtx, axis=1)),
            index=self.otus,
            columns=["prob"]
        )
        self.ratios_are_set = True
        if self.debug:
            print("Ratios are set")
    
    def noiser(self, seq, qual):
        bases = list("ATCGURYKMSWBDHVN")
        seq = list(seq)
        for i in range(len(seq)):
            prob = 10 ** (-qual[i] / 10) # 10 ^ (-Q / 10)
            rplc = [seq[i]] + [base for base in bases if base != seq[i]]
            probs = [1-prob] + [prob / len(rplc) for _ in range(len(rplc))]
            new_base = np.random.choice(rplc, p=probs)
            seq[i] = new_base if new_base != seq[i] else seq[i]            
        return "".join(seq)
    
    def generate(
            self,
            filename,
            reads=100,
            phred_params=(34, 5),
            ratios=None
        ):
        if self.debug:
            print(f"Generating {filename} sample")
        if not self.ratios_are_set:
            if ratios is None:
                raise ValueError(
                    "Ratios not set. Please set ratios before generating."
                )
            self.set_ratios(ratios)
        if ratios is not None:
            self.set_ratios(ratios)
        if not filename.endswith(".fastq"):
            raise ValueError("Output filename must be .fastq")
        chosen_otus = np.random.choice(
            self.otus,
            size=(reads,),
            p=self.otu_prob["prob"].values
        )
        chosen_regions = np.array([
            np.random.choice(
                self.regions,
                size=(reads,),
                p=softmax(self.ratio_ws.loc[otu, :].values)
            )
            for otu in chosen_otus
        ])
        chosen_idxs = [
            self.metadata[i]["idx"]
            for i in range(len(self.metadata))
            if self.metadata[i]["otu"] in chosen_otus
            and self.metadata[i]["region"] in chosen_regions
        ]
        records = [
            record
            for i, record in enumerate(SeqIO.parse(self.fraglib, "fasta"))
            if i in chosen_idxs
        ]
        
        for record in records:
            seq = str(record.seq)
            qual = np.int64(np.floor(phred_params[1] * np.random.randn(
                len(seq)
            ) + phred_params[0]))
            seq = self.noiser(seq, qual)
            record.seq = Seq(seq)
            record.letter_annotations = dict(
                phred_quality=qual
            )
        with open(filename, "w") as out:
            SeqIO.write(records, out, "fastq")
        print(f"File {filename} generated.")
        


def script(n_files=30, reads=500, path="./work", debug=False):
    import json, os
    path = os.path.abspath(path)
    if os.path.exists(path):
        print(f"Path {path} already exists, overwriting...")
    else:
        os.mkdir(path)
    fragment_library = "fragments.fasta"
    generator = FragmentGenerator(fragment_library, debug=debug)
    metadata = generator.metadata
    len_library = generator.fraglib_len
    with open(f"{path}/metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    otus = generator.otus
    file_mat = [
        softmax(np.random.randn(len(otus)))
        for _ in range(n_files)
    ]
    for i in range(n_files):
        fullpath = os.path.join(path, f"sample_{i}.fastq")
        print(f"Generating file {i+1}/{n_files}: {fullpath}")
        generator.generate(
            fullpath,
            reads=reads,
            ratios=file_mat[i]
        )
    file_mat = pd.DataFrame(
        file_mat,
        index=[f"sample_{i}" for i in range(n_files)],
        columns=otus
    )
    file_mat.to_csv(
        os.path.join(path, "random_otus_ratios.csv"),
        sep=",",
        header=True,
        index=True
    )
    print("Random reads generated")


if __name__ == "__main__":
    script(debug=True)