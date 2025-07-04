from generate_database import (
    align_matcher,
    fuzzy_matcher,
    exact_matcher,
    find_fragment,
    MATCHERS,
    LOG_LEVELS,
)
import typing as tp
from logger import (
    setup_logging,
    log_step,
)
from utils import (
    hash_sequence,
    build_primer_dict,
)
from concurrent.futures import ProcessPoolExecutor, as_completed
from tabulate import tabulate
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw
import numpy as np
import itertools
import joblib
import tempfile
import logging
import json
import os


def reducer_worker(i, j, seq_i, seq_j):
    alignment, score, _ = local_pairwise_align_ssw(
        DNA(str(seq_i.seq)),
        DNA(str(seq_j.seq)),
    )
    dist = 1 - (score / (2 * alignment.shape.position)) \
        if alignment.shape.position > 0 else 1.0
    return i, j, dist


def get_centroid_sequence(
    fasta_file: str,
    n_jobs: int = 1,
) -> SeqRecord:
    fasta_seqs = list(SeqIO.parse(fasta_file, "fasta"))
    n_seqs = len(fasta_seqs)
    combos = list(itertools.combinations(range(n_seqs), 2))
    dists = np.zeros((n_seqs, n_seqs), dtype=np.float32)
    if n_jobs > 1:
        with joblib.Parallel(n_jobs=n_jobs) as parallel:
            results = parallel(
                joblib.delayed(reducer_worker)(
                    i, j, fasta_seqs[i], fasta_seqs[j]
                ) for i, j in combos
            )
        for i, j, dist in results:
            dists[i, j] = dist
            dists[j, i] = dist
    else:
        for i, j in combos:
            dists[i, j] = reducer_worker(i, j, fasta_seqs[i], fasta_seqs[j])[2]
            dists[j, i] = dists[i, j]
    centroid_index = np.argmin(np.mean(dists, axis=1))
    return fasta_seqs[centroid_index]


@log_step
def build_centroid_database(
    fasta_dir: str,
    n_jobs: int = 1,
) -> str:
    fasta_dir = os.path.abspath(fasta_dir)
    fasta_seqs_files = [
        os.path.join(fasta_dir, f)
        for f in os.listdir(fasta_dir)
        if f.endswith(".fasta") and '16S' in f
    ]
    tmpfile = tempfile.mkstemp(suffix=".fasta")[1]
    writer = open(tmpfile, "a")
    for fasta_seq_file in fasta_seqs_files:
        sequence = get_centroid_sequence(
            fasta_file=fasta_seq_file,
            n_jobs=n_jobs
        )
        metadata = dict(
            file=os.path.splitext(
                os.path.basename(fasta_seq_file)
            )[0],
            seqhash=hash_sequence(str(sequence.seq)),
            seqid=sequence.id,
            seqname=sequence.name,
            seqdesc=sequence.description,
        )
        metadata["strain"] = "_".join(metadata["file"].split("_")[:2])
        metadata["name"] = "centroid_{strain}_{hash}_{seqid}".format(
            strain=metadata["strain"],
            hash=metadata["seqhash"],
            seqid=metadata["seqid"]
        )
        sequence.id = metadata["seqhash"]
        sequence.name = metadata["name"]
        sequence.description = json.dumps(metadata)
        SeqIO.write(sequence, writer, "fasta")
    writer.close()
    return tmpfile


def fragment_worker(
    fasta_record: SeqRecord,
    primer_dict: tp.Dict[str, Seq],
    region: str,
    matcher: str = "align_matcher",
    threshold: float = 0.8,
    reverse_primer3: bool = True,
):
    primer5 = primer_dict.get("F")
    primer3 = primer_dict.get("R")
    metadata = json.loads(fasta_record.description[21:])
    metadata["region"] = region
    fragment = find_fragment(
        sequence=fasta_record,
        primer5=primer5,
        primer3=primer3,
        matcher=MATCHERS[matcher],
        threshold=threshold,
        reverse_primer3=reverse_primer3
    )
    if isinstance(fragment, type(None)):
        return None
    return SeqRecord(
        fragment,
        id=metadata["seqhash"],
        name=metadata["name"],
        description=json.dumps(metadata)
    )


@log_step
def get_fragments_zymo(
    fasta_file: str,
    primer_file: str,
    outfile: str,
    matcher: str = "align_matcher",
    threshold: float = 0.8,
    reverse_primer3: bool = True,
    n_jobs: int = 1,
):
    primer_dict = build_primer_dict(primer_file)
    fasta_sequences = SeqIO.parse(fasta_file, "fasta")
    regions = list(primer_dict.keys())
    outfile = outfile if os.path.splitext(outfile)[1] == ".fasta" \
        else f"{outfile}.fasta"
    outfile = os.path.abspath(outfile)
    writer = open(outfile, "a")
    if n_jobs == 1:
        for fasta_sequence, region in itertools.product(fasta_sequences, regions):
            fragment = fragment_worker(
                fasta_record=fasta_sequence,
                primer_dict=primer_dict[region],
                matcher=matcher,
                threshold=threshold,
                reverse_primer3=reverse_primer3
            )
            if fragment is not None:
                SeqIO.write(fragment, writer, "fasta")
    else:
        kws_list = []
        for fasta_sequence, region in itertools.product(fasta_sequences, regions):
            kws_list.append(dict(
                fasta_record=fasta_sequence,
                primer_dict=primer_dict[region],
                region=region,
                matcher=matcher,
                threshold=threshold,
                reverse_primer3=reverse_primer3
            ))
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(
                    fragment_worker,
                    **kws
                ): kws for kws in kws_list
            }
            for future in as_completed(futures):
                fragment = future.result()
                if fragment is not None:
                    SeqIO.write(fragment, writer, "fasta")
    writer.close()


def main():
    import shutil
    setup_logging(log_level=LOG_LEVELS["DEBUG"])
    fasta_dir = "/data/16SDb/ZymoBIOMICS.STD.refseq.v2"
    logging.info(f"Building centroid database from {fasta_dir}")
    reduced_db = build_centroid_database(
        fasta_dir=os.path.join(fasta_dir, "ssrRNAs"),
        n_jobs=10,
    )
    shutil.copyfile(reduced_db, os.path.join(fasta_dir, "centroid.fasta"))
    outfile = os.path.join(fasta_dir, "fragments.fasta")
    logging.info(f"Saving fragment sequences to {outfile}")
    get_fragments_zymo(
        fasta_file=reduced_db,
        primer_file="/data/16SDb/primers.fasta",
        outfile=outfile,
        matcher="align_matcher",
        threshold=0.8,
        reverse_primer3=True,
        n_jobs=10,
    )

if __name__ == "__main__":
    main()
