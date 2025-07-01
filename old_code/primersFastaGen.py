from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

PRIMERS = {
    "V1": [
        "AGAGTTTGATCMTGGCTCAG",
        "TTACTCACCCGTICGCCRCT"
    ],
    "V2": [
        "AGYGGCGIACGGGTGAGTAA",
        "CYIACTGCTGCCTCCCGTAG"
    ],
    "V3": [
        "ACTCCTACGGGAGGCAGCAG",
        "ATTACCGCGGCTGCTGG"
    ],
    "V4": [
        "TGCCAGCAGCCGCGGTAA",
        "GGACTACARGGTATCTAAT"
    ],
    "V5": [
        "ATTAGATACCYTGTAGTCC",
        "CCGTCAATTCMTTTGAGTTT"
    ],
    "V6": [
        "AAACTCAAAKGAATTGACGG",
        "ACGAGCTGACGACARCCATG"
    ],
    "V7/V8": [
        "CATGGYTGTCGTCAGCTCGT",
        "ACGGGCGGTGTGTAC"
    ],
    "V9": [
        "GTACACACCGCCCGT",
        "TACCTTGTTACGACTT"
    ]
}


def main(primers):
    records = []
    for region, (fwd, rev) in primers.items():
        fwd = Seq(fwd)
        rev = Seq(rev)
        records.append(
            SeqRecord(
                fwd,
                id=f"{region}_F",
                description=f"{region} forward primer"
            )
        )
        records.append(
            SeqRecord(
                rev,
                id=f"{region}_R",
                description=f"{region} reverse primer"
            )
        )
    SeqIO.write(records, "primers.fasta", "fasta")


if __name__ == '__main__':
    main(primers=PRIMERS)