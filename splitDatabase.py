from Bio import SeqIO
import json
import numpy as np
import pandas as pd
from itertools import tee
from FragmentGenerator import extract_metadata
import os
import joblib

def fetch_fragment(idx, fr, region_indices):
    return fr if idx in region_indices else None


def split(n_jobs: int = 1):
    databases = os.path.abspath("databases")
    fragments = list(SeqIO.parse("fragments.fasta", "fasta"))
    metadata = json.load(open("fragment_metadata.json"))
    regions = np.unique([meta["region"] for meta in metadata])
    otus = np.unique([meta["otu"] for meta in metadata])
    regions = {
        region: [
            idx for idx, meta in enumerate(metadata)
            if meta["region"] == region
        ]
        for region in regions
    }
    os.makedirs(databases, exist_ok=True)
    for region in regions.keys():
        print(f"Splitting {region}...")
        region_indices = regions[region]
        completeness = len(region_indices) / len(otus) * 100
        print(f"Region completeness {completeness:.2f}% ({len(region_indices)}/{len(otus)})")
        print("Fetching fragments...")
        frags = list(joblib.Parallel(n_jobs=n_jobs, verbose=1)(
            joblib.delayed(fetch_fragment)(idx, fr, region_indices)
            for idx, fr in enumerate(fragments)
        ))
        frags = [f for f in frags if f is not None]
        print("Extracting metadata...")
        frag_meta = list(joblib.Parallel(n_jobs=n_jobs, verbose=1)(
            joblib.delayed(extract_metadata)(f)
            for f in frags
        ))
        with open(os.path.join(databases, f"fragments_{region}.fasta"), "w") as f:
            SeqIO.write(frags, f, "fasta")
        with open(os.path.join(databases, f"fragment_metadata_{region}.json"), "w") as f:
            json.dump(frag_meta, f, indent=4)
    print("Splitting completed.")

    


if __name__ == "__main__":
    split(6)