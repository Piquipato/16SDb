from concurrent.futures import ProcessPoolExecutor, as_completed
from skbio.alignment import local_pairwise_align_ssw
from tabulate import tabulate
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from skbio import DNA
import typing as tp
import pandas as pd
import joblib as jb
import numpy as np
import sys, os, re
import itertools
import functools
import traceback
import tempfile
import builtins
import hashlib
import logging
import inspect
import pickle
import shutil
import json


def track_open(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        original_open = builtins.open
        def open_logger(file, *a, **k):
            opened = os.path.abspath(file)
            logging.debug(f"Opened file {opened}")
            return original_open(file, *a, **k)
        builtins.open = open_logger
        try:
            result = func(*args, **kwargs)
        finally:
            builtins.open = original_open
        return result
    return wrapper
    

def log_step(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        runner = track_open(func)
        logging.debug(f"Running step {func.__name__}")
        prototype = inspect.signature(func)
        logging.debug(f"Prototype: {func.__name__}{prototype}")
        bound = prototype.bind(*args, **kwargs)
        bound.apply_defaults()
        params = [
            (k, v) for k, v in bound.arguments.items()
            if k not in ["args", "kwargs"]
        ]
        ag = bound.arguments.get("args")
        kw = bound.arguments.get("kwargs")
        if ag:
            params += [
                (f"args_{i}", ag[i]) for i in range(len(ag))
            ]
        if kw:
            params += [(k, v) for k, v in kw.items()]
        table = tabulate(
            params, headers=["Parameter", "Value"], tablefmt="pretty"
        )
        logging.debug(f"Arguments:\n{table}")
        try:
            return runner(*args, **kwargs)
        except Exception as e:
            logging.error(f"Error while running step {func.__name__}")
            logging.error(f"Exception: {e}")
            logging.error(f"Traceback:\n{traceback.format_exc()}")
            raise e
    return wrapper
        

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
        r"(\d+) ([\w.]+) ([\w. ]+) (\w+;.+;) ([\w]+)",
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


def reducer_worker(i, j, seq_i, seq_j):
    alignment, score, _ = local_pairwise_align_ssw(
        DNA(str(seq_i.seq)),
        DNA(str(seq_j.seq)),
    )
    dist = 1 - (score / (2 * alignment.shape.position)) \
        if alignment.shape.position > 0 else 1.0
    return i, j, dist

@log_step
def build_reduced_database(
    fasta_file: str,
    n_jobs: int = 1,
):
    logging.info("Reducing database to centroid sequences for each taxonomy")
    fasta_seqs = SeqIO.parse(fasta_file, "fasta")
    cp0, cp1, cp2, cp3, cp4, fasta_seqs = itertools.tee(fasta_seqs, 6)
    n_seqs = sum(1 for _ in cp0)
    logging.info("Extracting sequences metadata")
    metadata = [extract_metadata(seq) for seq in cp1]
    taxons = np.array([meta.get("taxon", "undefined") for meta in metadata])
    utaxons = np.unique(taxons)
    logging.info("Grouping sequences by taxonomy")
    taxons = {
        taxon: np.where(taxons == taxon)[0].tolist()
        for taxon in utaxons if taxon != "undefined"
    }
    tax_combos = {
        taxon: list(itertools.combinations(taxons[taxon], 2))
        for taxon in taxons.keys()
    }
    logging.info("Calculating sequence similarities for each taxonomy")
    combos = [tup for sublist in tax_combos.values() for tup in sublist]
    combos = sorted(combos, key=lambda x: (x[0], x[1]))
    args = []
    for i, seq_i in enumerate(cp2):
        for j, seq_j in enumerate(cp3):
            if (i, j) in combos:
                args.append((i, j, seq_i, seq_j))
    max_calcs = n_seqs * (n_seqs + 1) / 2
    efficiency = (1 - len(combos) / max_calcs) * 100
    logging.info(f"Doing {len(combos)} calculations, with efficiency {efficiency}")
    logging.info(f"Using {n_jobs} cores")
    if n_jobs == 1:
        results = [reducer_worker(*arg) for arg in args]
    else:
        results = jb.Parallel(n_jobs=n_jobs, verbose='1')(
            jb.delayed(reducer_worker)(*arg) for arg in args
        )
    distances = np.zeros((n_seqs, n_seqs), dtype=np.float32)
    for i, j, dist in results:
        distances[i, j] = dist
        distances[j, i] = dist
    distances = pd.DataFrame(distances)
    logging.info("Finding centroid sequences for each taxonomy")
    taxon_avg_dists = {
        taxon: dict(zip(
            idx_list,
            np.mean(distances.iloc[idx_list, idx_list].to_numpy(), axis=1)
        ))
        for taxon, idx_list in taxons.items()
    }
    taxon_min_samples = {
        taxon: list(dists.keys())[np.argmin(list(dists.values()))]
        for taxon, dists in taxon_avg_dists.items()
    }
    to_store = sorted(list(taxon_min_samples.values()))
    _, temp_fasta_file = tempfile.mkstemp(suffix=".fasta")
    logging.info(f"Storing reduced database to temporary file {temp_fasta_file}")
    with open(temp_fasta_file, "a") as f:
        for i, seq in enumerate(cp4):
            if not i in to_store:
                continue
            SeqIO.write(seq, f, "fasta")
    return temp_fasta_file
        

def convert_ambiguous_to_regex(primer: str):
    REGEXER = json.load(open("pattern.json", "r"))
    regex = ""
    for i in range(len(primer)):
        regex += REGEXER[primer[i]]
    return regex


def exact_matcher(
    sequence: SeqRecord,
    primer: tp.Union[Seq, SeqRecord],
    threshold: float = 0.8,
    reverse: bool = False
) -> tp.Tuple[int, float]:
    if reverse:
        primer = primer.reverse_complement()
    primer_regex = re.compile(convert_ambiguous_to_regex(primer))
    match = primer_regex.search(str(sequence))
    if match and reverse:
        return match.start(), 1.0
    elif match:
        return match.end(), 1.0
    return None, None


def fuzzy_matcher(
    sequence: SeqRecord,
    primer: tp.Union[Seq, SeqRecord],
    threshold: float = 0.8,
    reverse: bool = False
) -> tp.Tuple[int, float]:
    if reverse:
        primer = primer.reverse_complement()
    REGEXER = json.load(open("pattern.json", "r"))
    primer_patterns = [re.compile(REGEXER[base]) for base in primer]
    best = None
    best_score = 0
    seq = str(sequence)
    plen = len(primer)
    for i in range(len(seq) - plen + 1):
        window = seq[i:i+plen]
        matches = sum(
            pat.fullmatch(base) is not None 
            for pat, base in zip(primer_patterns, window)
        )
        score = matches / plen
        if score > best_score and score >= threshold:
            best = (i, i+plen, score)
            best_score = score
    final_score = best[2] if best else 0
    if final_score < threshold:
        return None, None
    if reverse:
        return best[0], best[2]
    else:
        return best[1], best[2]


def align_matcher(
    sequence: SeqRecord,
    primer: tp.Union[Seq, SeqRecord],
    threshold: float = 0.8,
    reverse: bool = False
) -> tp.Tuple[int, float]:
    if reverse:
        primer = primer.reverse_complement()
    best_score = None
    best_aln = None
    for one_primer in sorted([
        str(s) for s in DNA(str(primer.seq)).expand_degenerates()
    ]):
        this_aln = local_pairwise_align_ssw(
            DNA(one_primer),
            DNA(str(sequence.seq)),
        )
        score = this_aln[1]
        if best_score is None or score > best_score:
            best_score = score
            best_aln = this_aln
    if best_aln is None:
        return None, None
    (aln_prim, aln_seq), score, (prim_pos, seq_pos) = best_aln
    match_start = max(0, seq_pos[0] - prim_pos[0])
    match_end = min(
        len(str(sequence.seq)),
        seq_pos[1] + len(primer) - prim_pos[1]
    )

    primer = str(primer.seq)
    bits_prim = [primer[:prim_pos[0]], aln_prim, primer[prim_pos[1]+1:]]
    full_prim = ''.join(map(str, bits_prim))
    sequence = str(sequence.seq)
    bits_seq = [
        '-' * match_start,
        sequence[match_start:seq_pos[0]],
        aln_seq,
        sequence[seq_pos[1]:match_end],
        '-'*(len(sequence)-match_end)
    ]
    full_seq = ''.join(map(str, bits_seq))

    matches = sum(
        s in DNA.degenerate_map.get(p, {p})
        for p, s in zip(full_prim, full_seq[match_start:match_end])
    )
    
    if reverse:
        amplicon_pos = match_start
    else:
        amplicon_pos = match_end
    final_score = matches / len(primer)
    if final_score < threshold:
        return None, None
    return amplicon_pos, final_score


def find_fragment(
    sequence: SeqRecord,
    primer5: tp.Union[Seq, SeqRecord],
    primer3: tp.Union[Seq, SeqRecord],
    matcher: tp.Callable = align_matcher,
    threshold: float = 0.8,
    reverse_primer3: bool = True
):
    if not reverse_primer3:
        primer3 = primer3.reverse_complement()
    seq = str(sequence.seq)
    match5, score5 = matcher(
        sequence,
        primer5,
        threshold=threshold,
        reverse=False
    )
    if match5 is None:
        return None
    match3, score3 = matcher(
        sequence,
        primer3,
        threshold=threshold,
        reverse=True
    )
    if match3 is None:
        return None
    return Seq(seq[match5:match3])


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
            primer_dict[region][orientation] = Seq(str(primer.seq))
    return {
        region: primers for region, primers in primer_dict.items()
        if "F" in primers.keys() and "R" in primers.keys()
    }


def fragment_worker(
    fasta_record: SeqRecord,
    primer_dict: tp.Dict[str, Seq],
    matcher: tp.Callable = align_matcher,
    threshold: float = 0.8,
    reverse_primer3: bool = True,
):
    primer5 = primer_dict.get("F")
    primer3 = primer_dict.get("R")
    metadata = extract_metadata(fasta_record)
    fragment = find_fragment(
        sequence=fasta_record,
        primer5=primer5,
        primer3=primer3,
        matcher=matcher,
        threshold=threshold,
        reverse_primer3=reverse_primer3
    )
    if isinstance(fragment, type(None)):
        return None
    return SeqRecord(
        fragment,
        id=metadata["seqhash"],
        name=metadata["uname"],
        description=json.dumps(metadata)
    )

@log_step
def build_fragment_database(
    fasta_file: str,
    primer_file: str,
    outfile: str,
    matcher: tp.Callable = align_matcher,
    threshold: float = 0.8,
    reverse_primer3: bool = True,
    n_jobs: int = 1
):
    logging.info("Starting to build fragment database")
    logging.info("Building region dictionary with primers...")
    primer_dict = build_primer_dict(primer_file)
    fasta_sequences = SeqIO.parse(fasta_file, "fasta")
    regions = list(primer_dict.keys())
    outfile = outfile if os.path.splitext(outfile)[0] == ".fasta" \
        else f"{outfile}.fasta"
    outfile = os.path.abspath(outfile)
    logging.info(f"Writing fragment database to file {outfile}")
    writer = open(outfile, "a")
    logging.info(f"Calculating fragments with {n_jobs} cores")
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
            sig = inspect.signature(fragment_worker)
            for future in as_completed(futures):
                kws = futures[future]
                try:
                    fragment = future.result()
                    if fragment is not None:
                        SeqIO.write(fragment, writer, "fasta")
                except Exception as e:
                    logging.error(f"Error while calculating fragment")
                    table = tabulate(
                        [(k, v) for k, v in kws.items()],
                        headers=["Parameter", "Value"]
                    )
                    logging.error(f"Function: {fragment_worker.__name__}{sig}")
                    logging.error(f"Arguments:\n{table}")
                    logging.error(f"Exception: {e}")
                    logging.error(f"Traceback:\n{traceback.format_exc()}")
                    continue
    writer.close()


def setup_logging(
    log_level: int = 20,
    log_file: str = None,
    log_format: str = '[%(asctime)s | %(levelname)s] : %(message)s'
):
    logger = logging.getLogger()
    logger.setLevel(log_level)
    formatter = logging.Formatter(log_format)
    if not any(
        isinstance(handler, logging.StreamHandler)
        for handler in logger.handlers
    ):
        stream = logging.StreamHandler()
        stream.setLevel(log_level)
        stream.setFormatter(formatter)
        logger.addHandler(stream)
    if not any(
        isinstance(handler, logging.StreamHandler)
        for handler in logger.handlers
    ) and not isinstance(log_file, type(None)):
        fileh = logging.FileHandler(filename=log_file)
        fileh.setLevel(log_level)
        fileh.setFormatter(formatter)
        logger.addHandler(fileh)


def main(
    fasta_file: str,
    primer_file: str,
    outfile: str,
    matcher: tp.Callable = align_matcher,
    threshold: float = 0.8,
    reverse_primer3: bool = True,
    n_jobs: int = 1,
    log_level: int = 20
):
    setup_logging(log_level=log_level)
    reduced_fasta = build_reduced_database(
        fasta_file=fasta_file,
        n_jobs=n_jobs,
    )
    build_fragment_database(
        fasta_file=reduced_fasta,
        primer_file=primer_file,
        outfile=outfile,
        matcher=matcher,
        threshold=threshold,
        reverse_primer3=reverse_primer3,
        n_jobs=n_jobs,
    )


zipper = lambda x: dict(zip([func.__name__ for func in x], x))
MATCHERS = zipper([
    align_matcher,
    exact_matcher,
    fuzzy_matcher,
])
LOG_LEVELS = dict(sorted([
    (name, lvl) for name, lvl in logging._nameToLevel.items()
], key=lambda x: x[1]))

def cli(*args):
    import argparse
    parser = argparse.ArgumentParser(
        prog="generate_database",
        description="Script used to generate fragment databases given " \
            "a set of forward (F) and reverse (R) primers and a file with " \
            "sequences."
    )
    parser.add_argument(
        "--fasta-file", "--file", "-f",
        dest="fasta_file",
        type=str,
        required=True,
        help="File in .fasta format with the sequences " \
            "to grab the fragments from."
    )
    parser.add_argument(
        "--primer-file", "--p-file", "-p",
        dest="primer_file",
        type=str,
        required=True,
        help="File with the primers to be used to extract " \
            "the fragments. Read docs on how to format this file (.fasta)."
    )
    default_output = os.path.abspath("./output.fasta")
    parser.add_argument(
        "--outfile", "-o", "--output-file",
        dest="outfile",
        type=str,
        default=default_output,
        help="File in .fasta format in which to store fragments " \
            f"(default: {default_output})."
    )
    parser.add_argument(
        "--matcher", "-m",
        dest="matcher",
        type=str,
        default="align_matcher",
        choices=list(MATCHERS.keys()),
        help="Matcher to use to search for the primers in the sequences " \
            "(default: align_matcher)."
    )
    parser.add_argument(
        "--threshold", "-t",
        dest="threshold",
        type=float,
        default=0.8,
        help="Threshold of similarity to keep matches of primers. " \
            "Bounded between 0 and 1."
    )
    parser.add_argument(
        "--no-reverse",
        dest="reverse_primer3",
        action="store_false",
        default=True,
        help="Include this flag if you do not wish to reverse-complement " \
            "the 3' primers sequences provided."
    )
    parser.add_argument(
        "--n-jobs", "--n-cpus", "-j",
        dest="n_jobs",
        type=int,
        default=1,
        help="Number of cores to use."
    )
    parser.add_argument(
        "--log-level",
        dest="log_level",
        type=str,
        default="INFO",
        choices=list(LOG_LEVELS.keys()),
        help="Log level to use (default: INFO)"
    )
    cliargs = parser.parse_args(*args)
    kwargs = dict(
        fasta_file=os.path.abspath(cliargs.fasta_file),
        primer_file=os.path.abspath(cliargs.primer_file),
        outfile=os.path.abspath(cliargs.outfile),
        matcher=MATCHERS[cliargs.matcher],
        threshold=cliargs.threshold,
        reverse_primer3=cliargs.reverse_primer3,
        n_jobs=cliargs.n_jobs,
        log_level=LOG_LEVELS[cliargs.log_level],
    )
    main(**kwargs)


if __name__ == "__main__":
    cli()