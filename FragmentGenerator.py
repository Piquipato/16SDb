from Bio.SeqRecord import SeqRecord
from scipy.special import softmax
from Bio.Seq import Seq
from Bio import SeqIO
import typing as tp
import pandas as pd
import joblib as jb
import numpy as np
import os, re, sys
import warnings
import json


"""
    TODO: Hacer Docker con entorno de ejecucion y dependencias
    TODO: Hacer un README
    TODO: Hacer un environment.yml
    TODO: Hacer un requirements.txt
    TODO: Ordenar repositorio
    x TODO: Hacer Forward y Reverse
    TODO: Dividir base de datos en fragmentos (V1-9)
    TODO: Mirar distribuciones y mirar cuaderno (fastqc y multiqc)
    TODO: Meter control de que region meter en los archivos
    x TODO: Meter parametros de longitud de fragmentos
    TODO: Medir distribuciones (calidad por base, calidad media fragmento, 
                                longitud, numero de lecturas)
    - TODO: Meter calidad por fragmentos (calidad media o calidad por base, flag)
    TODO: Mirar OTUs Angelo
"""


def extract_metadata(fragment: SeqRecord, mode: str = "normal") -> tp.Dict[str, str]:
    props = {"id": "","name": "","description": ""}
    for key in props.keys():
        props[key] = getattr(fragment, key)
    if mode == "json":
        dsc_json = json.loads(props["description"])
        props["region"] = dsc_json.get("region", "")
        props["ggid"] = dsc_json.get("id", "")
        props["ggname"] = dsc_json.get("name", "")
        props["ggdesc"] = dsc_json.get("description", "")
    else:
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
            r"(\w{1}__[\w]+;|Unclassified;)",
            props["ggdesc"]
        ) if m is not None
    ])
    props["ggtag"] = props["ggdesc"].split(" ")[1]
    return props

def extract_metadata_library(library, mode = "normal"):
    metadata = []
    len_library = len(library)
    for i, f in enumerate(SeqIO.parse(library, "fasta")):
        props = extract_metadata(f, mode=mode)
        props["idx"] = i
        metadata.append(props)
    return metadata, len_library

class FragmentGenerator:

    def __init__(
        self,
        fasta_path: str,
        metadata_path: str = None,
        debug: bool = False,
    ):
        self.debug = debug
        if self.debug:
            print("Initializing FragmentGenerator...")
        self.fasta_path = fasta_path
        if metadata_path is None:
            if self.debug:
                print("Extracting metadata from fasta library...")
            self.metadata, _ = extract_metadata_library(fasta_path, mode="normal")
        else:
            if self.debug:
                print("Loading metadata from json file...")
            self.metadata = json.load(open(metadata_path))
        self.otus_idx = [props["otu"] for props in self.metadata]
        self.regions_idx = [props["region"] for props in self.metadata]
        self.otus = np.unique(self.otus_idx)
        self.regions = np.unique(self.regions_idx)
        self.present_matrix = pd.DataFrame(np.zeros(
            shape=(len(self.otus), len(self.regions)),
            dtype=np.int64
        ), index=self.otus, columns=self.regions)
        completeness = np.sum(self.present_matrix.to_numpy().flatten()) / \
            np.size(self.present_matrix.to_numpy().flatten()) * 100
        if self.debug:
            print("Completeness percentage of database: " \
                  f"{completeness:.2f}%")
            print("FragmentGenerator initialized.")
    
    def noiser(
        self,
        record: SeqRecord,
        phred_params: tp.Tuple[int, int] = (34, 5),
        length_params: tp.Tuple[int, int] = (200, 30),
        reverse_complement: bool = False
    ):
        length = np.min([
            np.int64(np.floor(
                length_params[1] * np.random.randn() + length_params[0]
            )),
            len(record.seq)
        ])
        idx = np.random.randint(0, len(record.seq) - length) \
            if length != len(record.seq) else 0
        if reverse_complement:
            seq = list(str(record.seq.reverse_complement()))[idx:idx+length]
        else:
            seq = list(str(record.seq))[idx:idx+length]
        qual = np.int64(np.floor(phred_params[1] * np.random.randn(
            len(seq)
        ) + phred_params[0]))
        bases = list("ATCGURYKMSWBDHVN")
        for i in range(len(seq)):
            prob = 10 ** (-qual[i] / 10) # 10 ^ (-Q / 10)
            rplc = [seq[i]] + [base for base in bases if base != seq[i]]
            probs = [1-prob] + [prob / (len(rplc)-1) for _ in range(len(rplc)-1)]
            new_base = np.random.choice(rplc, p=probs)
            seq[i] = new_base if new_base != seq[i] else seq[i]
        record.seq = Seq("".join(seq))
        record.letter_annotations = dict(phred_quality=qual)
        return record

    def choose_fragments(
            self,
            n_reads: int = 100,
            ratios: tp.Iterable[float] = None,
            phred_params: tp.Tuple[int, int] = (34, 5),
            length_params: tp.Tuple[int, int] = (200, 30),
            reverse_complement: bool = False
        ):
        chosen_otus = np.random.choice(
            self.otus,
            size=(n_reads,),
            p=ratios
        )
        chosen_regions = [
            np.random.choice(
                np.array(self.present_matrix.columns)[np.where(
                    self.present_matrix.loc[chosen_otus[i], :] == 1
                )[0]]
            )
            for i in range(len(chosen_otus))
        ]
        chosen_idxs = [
            i for i, meta in enumerate(self.metadata)
            if meta["otu"] in chosen_otus
            and meta["region"] in chosen_regions
        ]
        if self.debug:
            print("Grabbing records and noising them...")
        records = [
            self.noiser(
                record,
                phred_params=phred_params,
                length_params=length_params,
                reverse_complement=reverse_complement
            )
            for i, record in enumerate(SeqIO.parse(self.fasta_path, "fasta"))
            if i in chosen_idxs
        ]
        return records


    def generate(
            self,
            filename: str,
            ratios: tp.Iterable[float] = None,
            n_reads: int = 100,
            phred_params: tp.Tuple[int, int] = (34, 5),
            length_params: tp.Tuple[int, int] = (200, 30),
            forward_reverse: bool = False,
        ):
        if ratios is None:
            if self.debug:
                print("No ratios provided, generating random ratios...")
            ratios = softmax(np.random.randn(len(self.otus)))
        else:
            if self.debug:
                print("Using provided ratios...")
            ratios = np.array(ratios)
            if len(ratios) != len(self.otus):
                raise ValueError("Ratios must be the same length as OTUs")
            if np.sum(ratios) != 1:
                ratios = softmax(ratios)
        if not filename.endswith(".fastq"):
            raise ValueError("Output filename must be .fastq")
        if self.debug:
            print("Choosing random fragments...")
        if forward_reverse:
            name, ext = os.path.splitext(filename)
            filenames = [
                (f"{name}_R1{ext}", False),
                (f"{name}_R2{ext}", True)
            ]
        else:
            filenames = [(filename, False)]
        for file, reverse_complement in filenames:
            if self.debug:
                print(f"Generating fragments for sample {file}...")
            records = self.choose_fragments(
                n_reads=n_reads,
                ratios=ratios,
                phred_params=phred_params,
                length_params=length_params,
                reverse_complement=reverse_complement
            )
            if self.debug:
                print(f"Writing records to file {file}...")
            with open(file, "w") as out:
                SeqIO.write(records, out, "fastq")
            if self.debug:
                print(f"Sample {file} generated.")
        return ratios


def script(
    fasta_path: str,
    output_dir: str = "./work",
    metadata_path: str = None,
    ratios_mtx_path: str = None,
    get_otus_list: bool = False,
    n_reads: int = 100,
    n_samples: int = 30,
    phred_params: tp.Tuple[int, int] = (34, 5),
    length_params: tp.Tuple[int, int] = (200, 30),
    forward_reverse: bool = False,
    n_jobs: int = 1,
    debug: bool = False,
):
    if debug:
        print("Starting script to create synthetic reads...")
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    elif os.path.isfile(output_dir):
        raise ValueError(f"Output directory {output_dir} is a file.")
    else:
        warnings.warn(f"Output directory {output_dir} already exists." \
                      " Overwriting files...")
    if not os.path.isfile(fasta_path) or not fasta_path.endswith(".fasta"):
        raise ValueError("Input file must be a fasta file.")

    fg = FragmentGenerator(
        fasta_path,
        metadata_path=metadata_path,
        debug=debug
    )
    if get_otus_list:
        pd.DataFrame(
            np.zeros_like(fg.otus),
            columns=["sample_1"],
            index=list(fg.otus)
        ).to_csv(
            os.path.join(output_dir, "otus.csv"),
            index=True,
            header=True
        )
        return 0
    if ratios_mtx_path is None:
        name_files = [f"sample_{i + 1}.fastq" for i in range(n_samples)]
        ratios_mtx = dict(zip(
            name_files,
            [None for _ in range(n_samples)]
        ))
    else:
        if not os.path.isfile(ratios_mtx_path) or \
            not ratios_mtx_path.endswith(".csv"):
            raise ValueError("Ratios matrix must be a csv file.")
        ratios_mtx = pd.read_csv(
            ratios_mtx_path,
            sep=",",
            index_col=0,
            header=0
        )
        if ratios_mtx.shape[0] != len(fg.otus):
            raise ValueError("Number of OTUs in matrix does not match.")
        name_files = list(ratios_mtx.columns)
        n_reads = ratios_mtx.shape[1]
    filenames = [
        os.path.join(output_dir, name)
        for name in name_files
    ]
    if debug:
        print("Generating samples...")
    if n_jobs == 1:
        ratios = []
        for filename in filenames:
            ratios.append(fg.generate(
                filename,
                ratios=ratios_mtx[os.path.basename(filename)],
                n_reads=n_reads,
                phred_params=phred_params,
                length_params=length_params,
                forward_reverse=forward_reverse
            ))
    else:
        ratios = jb.Parallel(n_jobs=n_jobs, verbose=debug)(
            jb.delayed(fg.generate)(
                filename,
                ratios=ratios_mtx[os.path.basename(filename)],
                n_reads=n_reads,
                phred_params=phred_params,
                length_params=length_params,
                forward_reverse=forward_reverse
            ) for filename in filenames
        )
    ratios = pd.DataFrame(
        np.array(ratios),
        columns=fg.otus,
        index=name_files
    ).transpose() 
    ratios.to_csv(
        os.path.join(output_dir, "ratios.csv"),
        index=True,
        header=True
    )
    return 0

def cli():
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate synthetic reads from a fasta file."
    )
    parser.add_argument(
        "--fasta-path", "--fasta", "-f",
        required=True,
        dest="fasta_path",
        type=str,
        help="Path to the fasta file."
    )
    parser.add_argument(
        "-o", "--output-dir",
        dest="output_dir",
        type=str,
        default="./work",
        help="Path to the output directory."
    )
    parser.add_argument(
        "-m", "--metadata-path", "--metadata",
        dest="metadata_path",
        type=str,
        default=None,
        help="Path to the metadata file."
    )
    parser.add_argument(
        "-r", "--ratios-mtx-path", "--ratios",
        dest="ratios_mtx_path",
        type=str,
        default=None,
        help="Path to the ratios matrix file."
    )
    parser.add_argument(
        "-n", "--n-reads",
        dest="n_reads",
        type=int,
        default=100,
        help="Number of reads to generate."
    )
    parser.add_argument(
        "-s", "--n-samples",
        dest="n_samples",
        type=int,
        default=30,
        help="Number of samples to generate."
    )
    parser.add_argument(
        "-nfr", "--no-forward-reverse",
        dest="forward_reverse",
        action="store_false",
        default=True,
        help="Generate forward and reverse reads."
    )
    parser.add_argument(
        "-p", "--phred-params",
        dest="phred_params",
        type=str,
        default="(34, 5)",
        help="Phred parameters for noise generation (Tuple[int, int])."
    )
    parser.add_argument(
        "-l", "--length-params",
        dest="length_params",
        type=str,
        default="(100, 15)",
        help="Length parameters for fragment generation (Tuple[int, int])."
    )
    parser.add_argument(
        "-j", "--n-jobs",
        dest="n_jobs",
        type=int,
        default=1,
        help="Number of jobs to run in parallel."
    )
    parser.add_argument(
        "-g", "--get-otus-list",
        dest="get_otus_list",
        action="store_true",
        default=False,
        help="Get the list of OTUs."
    )
    parser.add_argument(
        "-d", "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help="Enable debug mode."
    )
    args = parser.parse_args()
    script(
        fasta_path=args.fasta_path,
        output_dir=args.output_dir,
        metadata_path=args.metadata_path,
        ratios_mtx_path=args.ratios_mtx_path,
        get_otus_list=args.get_otus_list,
        n_reads=args.n_reads,
        n_samples=args.n_samples,
        phred_params=tuple(
            int(it) for it in args.phred_params.strip("()").split(",")
        ),
        length_params=tuple(
            int(it) for it in args.length_params.strip("()").split(",")
        ),
        forward_reverse=args.forward_reverse,
        n_jobs=args.n_jobs,
        debug=args.debug
    )


if __name__ == "__main__":
    cli()
    