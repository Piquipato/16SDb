from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import typing as tp
import hashlib
import pickle
import os, re
import json
import os


def hash_sequence(sequence, size=10):
    h = hashlib.blake2b(digest_size=size)
    h.update(pickle.dumps(sequence))
    return h.hexdigest()


def extract_metadata(sequence: SeqRecord) -> tp.Dict[str, str]:
    metadata = dict(
        seqid=sequence.id,
        seqname=sequence.name,
        seqhash=hash_sequence(sequence.seq),
        description=sequence.description,
    )
    metadata["uname"] = "sequence_{id}_{name}_{hash}".format(
        id=metadata["seqid"],
        name=metadata["seqname"],
        hash=metadata["seqhash"],
    )
    search = re.search(
        r"(\d+) ([\w.]+) ([\w. ]+) (Archaea;.+;|Bacteria;.+;) ([\w]+)",
        metadata.get("description", "")
    )
    for i, key in enumerate([
        "dscid",
        "strnid",
        "strnname",
        "taxon",
        "otu",
    ]):
        metadata[key] = search.group(i + 1) if search else "undefined"
    return metadata


def convert_ambiguous_to_regex(primer: str):
    REGEXER = json.load(open(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "data/pattern.json"))
        , "r"))
    regex = ""
    for i in range(len(primer)):
        regex += REGEXER[primer[i]]
    return regex


def build_primer_dict(
    primer_file: str,
) -> tp.Dict[str, tp.Dict[str, Seq]]:
    primer_sequences = SeqIO.parse(primer_file, "fasta")
    primer_dict = {}
    for primer in primer_sequences:
        region = "_".join(primer.name.split("_")[:-1])
        orientation = primer.name.split("_")[-1]
        if not region in primer_dict.keys():
            primer_dict[region] = {}
        if orientation not in primer_dict[region].keys():
            primer_dict[region][orientation] = primer.seq
    return {
        region: primers for region, primers in primer_dict.items()
        if "F" in primers.keys() and "R" in primers.keys()
    }

