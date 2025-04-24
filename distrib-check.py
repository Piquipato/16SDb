from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import typing as tp
import matplotlib.pyplot as plt
from itertools import tee, islice
from concurrent.futures import ProcessPoolExecutor, as_completed
import joblib as jb
import pandas as pd
import numpy as np
import os

def batched(it, n):
    while (batch := tuple(islice(it, n))):
        yield batch

class IncrementalHistogram:

    def __init__(
        self,
        data: tp.Union[tp.Iterable[float], tp.Iterable[int]],
        bins: int = 100,
        n_jobs: int = 1
    ):
        self.n_jobs = n_jobs
        cp1, cp2, self.data = tee(data, 3)
        self.range = self.find_range(cp1)
        self.discrete = self.find_type(cp2) == int
        self.bins = np.linspace(self.range[0], self.range[1], bins + 1) \
            if not self.discrete \
            else np.linspace(
                self.range[0], self.range[1],
                len(range(self.range[0], self.range[1] + 1))
            )
        self.histogram = np.zeros(len(self.bins) - 1, dtype=np.float64)
        if self.n_jobs == 1:
            for pt in self.data:
                self.histogram += np.histogram(
                    [pt], bins=self.bins, range=self.range
                )[0]
        else:
            with ProcessPoolExecutor(max_workers=self.n_jobs) as executor:
                futures = {
                    executor.submit(self.calc_histogram_update, pt): pt
                    for pt in self.data
                }
                for future in as_completed(futures):
                    try:
                        result = future.result()
                        self.histogram += result
                    except Exception as exc:
                        print(f"Error: {exc}")
                        continue
    
    def __call__(self) -> tp.Tuple[np.ndarray, np.ndarray]:
        return self.bins, self.histogram

    def calc_histogram_update(self, pt):
        return np.histogram([pt], bins=self.bins, range=self.range)[0]

    def find_xtr_chunk(
        self,
        chunk: tp.Union[tp.Iterable[float], tp.Iterable[int]],
        minmax: bool = False
    ):
        return max(chunk) if minmax else min(chunk)
    
    def parallel_xtr(
        self,
        array: tp.Union[tp.Iterable[float], tp.Iterable[int]],
        minmax: bool = False,
        n_jobs: int = 1
    ) -> float:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(self.find_xtr_chunk, chunk, minmax): chunk
                for chunk in batched(array, n_jobs)
            }
            results = []
            for future in as_completed(futures):
                try:
                    result = future.result()
                    results.append(result)
                except Exception as exc:
                    print(f"Error: {exc}")
                    continue
        return max(results) if minmax else min(results)

    def find_range(
        self,
        data: tp.Union[tp.Iterable[float], tp.Iterable[int]]
    ) -> tp.Tuple[float, float]:
        d1, d2 = tee(data)
        if self.n_jobs == 1:
            return min(d1), max(d2)
        else:
            min_val = self.parallel_xtr(array=d1, minmax=False, n_jobs=self.n_jobs)
            max_val = self.parallel_xtr(array=d2, minmax=True, n_jobs=self.n_jobs)
            return min_val, max_val
            
    
    def find_type(
        self,
        data: tp.Union[tp.Iterable[float], tp.Iterable[int]]
    ) -> tp.Type:
        cp1, cp2, cp3 = tee(data, 3)
        if any(
            jb.Parallel(n_jobs=self.n_jobs)
            (
                jb.delayed(
                    lambda i: isinstance(i, float) \
                    or isinstance(i, np.float32)
                    or isinstance(i, np.float64)
                )(i) for i in cp1
            )
        ):
            return float
        if all(
            jb.Parallel(n_jobs=self.n_jobs)
            (
                jb.delayed(
                    lambda i: isinstance(i, int) \
                    or isinstance(i, np.int32)
                    or isinstance(i, np.int64)
                )(i) for i in cp2
            )
        ):
            return int
        types = set(
            jb.Parallel(n_jobs=self.n_jobs)
            (
                jb.delayed(
                    lambda i: type(i)
                )(i) for i in cp3
            )
        ).difference(set([
            int,
            np.int32,
            np.int64,
            float,
            np.float32,
            np.float64
        ]))
        raise TypeError(f"Anomalous data types: {types}")



def build_samples_iterator(
    dirname: str,
) -> tp.Iterator[str]:
    for filename in os.listdir(dirname):
        if filename.endswith(".fastq"):
            filepath = os.path.join(dirname, filename)
            yield filepath


def build_record_iterator(
    dirname: str
) -> tp.Iterator[SeqRecord]:
    samples = build_samples_iterator(dirname)
    for filepath in samples:
        with open(filepath, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                yield record


def check_base_quality_distribution(
    records: tp.Iterable[SeqRecord],
    n_jobs: int = 1
) -> IncrementalHistogram:
    print("Counting base qualities...")
    def build_base_quality_iterator(
        records: tp.Iterable[SeqRecord],
    ) -> tp.Iterator[int]:
        for record in records:
            phred_scores = record.letter_annotations["phred_quality"]
            for score in phred_scores:
                yield score
    phred_scores = build_base_quality_iterator(records)
    return IncrementalHistogram(phred_scores, n_jobs=n_jobs)


def check_fragment_mean_quality_distribution(
    records: tp.Iterable[SeqRecord],
    n_jobs: int = 1
) -> IncrementalHistogram:
    print("Counting fragment mean qualities...")
    def build_fragment_mean_quality_iterator(
        records: tp.Iterable[SeqRecord],
    ) -> tp.Iterator[float]:
        for record in records:
            mean_quality = np.mean(record.letter_annotations["phred_quality"])
            yield mean_quality
    mean_phred_scores = build_fragment_mean_quality_iterator(records)
    return IncrementalHistogram(mean_phred_scores, n_jobs=n_jobs)


def check_fragment_length_distribution(
    records: tp.Iterable[SeqRecord],
    n_jobs: int = 1
) -> IncrementalHistogram:
    print("Counting fragment lengths...")
    def build_fragment_length_iterator(
        records: tp.Iterable[SeqRecord],
    ) -> tp.Iterator[int]:
        for record in records:
            yield len(record.seq)
    fragment_lengths = build_fragment_length_iterator(records)
    return IncrementalHistogram(fragment_lengths, n_jobs=n_jobs)


def check_read_number_distribution(
    samples: tp.Iterable[str],
    n_jobs: int = 1
) -> IncrementalHistogram:
    print("Counting reads...")
    def build_read_number_iterator(
        samples: tp.Iterable[str],
    ) -> tp.Iterator[int]:
        for sample in samples:
            read_count = 0
            for read in SeqIO.parse(sample, "fastq"):
                read_count += 1
            yield read_count
    read_counts = build_read_number_iterator(samples)
    return IncrementalHistogram(read_counts, n_jobs=n_jobs)


def get_distributions_from_database(
    db_path: str,
    n_jobs: int = 1
) -> tp.Dict[str, IncrementalHistogram]:
    print("Getting iterators...")
    samples = build_samples_iterator(db_path)
    r1, r2, r3 = tee(build_record_iterator(db_path), 3)
    return {
        "base_quality": check_base_quality_distribution(r1, n_jobs=n_jobs),
        "fragment_quality": check_fragment_mean_quality_distribution(r2, n_jobs=n_jobs),
        "fragment_length": check_fragment_length_distribution(r3, n_jobs=n_jobs),
        "read_number": check_read_number_distribution(samples, n_jobs=n_jobs)
    }


def check_distributions(n_jobs: int = 1) -> None:
    print("Checking distributions...")
    dbs_path = os.path.abspath("fastq")
    distrib_path = os.path.join(os.path.split(dbs_path)[0], "distributions")
    print(f"Distributions will be saved to {distrib_path}")
    os.makedirs(distrib_path, exist_ok=True)
    for db_name in os.listdir(dbs_path):
        print(f"Checking {db_name}...")
        os.makedirs(os.path.join(distrib_path, db_name), exist_ok=True)
        db_path = os.path.join(dbs_path, db_name)
        distributions = get_distributions_from_database(db_path, n_jobs=n_jobs)
        for name, distribution in distributions.items():
            print(f"Generating {name} distribution for {db_name} plot and data...")
            bins, histogram = distribution()
            plt.figure(figsize=(10, 6))
            plt.bar(bins[:-1], histogram, width=np.diff(bins), align='edge')
            plt.title(f"{name} distribution")
            plt.xlabel(name)
            plt.ylabel("frequency")
            plt.grid()
            plt.savefig(os.path.join(distrib_path, f"{name}_distribution.png"))
            plt.close()
            pd.DataFrame(
                np.transpose(np.array([bins[:-1], histogram])),
                columns=["bins", "frequency"]
            ).to_csv(
                os.path.join(distrib_path, f"{name}_distribution.csv")
            )



if __name__ == "__main__":
    check_distributions(n_jobs=5)

