from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

"""
TO-DOs:
o Escribir una función que convierta una secuencia de nucleótidos con 
  ambigüedades a una expresión regular.
o Iterar sobre todos los primers presentes en el archivo primers.fasta.
o Iterar sobre todas las secuencias presentes en el archivo.
o Paralelizar la iteración sobre las secuencias del archivo.
o Guardar los fragmentos en un archivo fasta nuevo.
o Crear un pd.DataFrame con (bicho, región 16S) == fragmento.
o Guardar el DataFrame en un archivo csv/xlsx.
"""


def convert_ambiguous_to_regex(primer):
    primer_seq = Seq(primer, IUPAC.ambiguous_dna)
    return primer_seq.tomutable().toseq().ungap("-").back_transcribe().translate(to_stop=True).tostring()


def find_fragment(sequence, primer5, primer3):
    primer5_regex = re.compile(convert_ambiguous_to_regex(primer5))
    primer3_regex = re.compile(convert_ambiguous_to_regex(primer3))

    match5 = primer5_regex.search(sequence)
    match3 = primer3_regex.search(sequence, match5.end() if match5 else 0)

    if match5 and match3:
        return sequence[match5.end():match3.start()]
    return None


def main(fasta_file, primer_file):
    # Load fasta file
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    primer_sequences = list(SeqIO.parse(open(primer_file), 'fasta'))

    primer5 = str(primer_sequences[0].seq)
    primer3 = str(primer_sequences[1].seq).reverse_complement()

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fragment = find_fragment(sequence, primer5, primer3)
        if fragment:
            print(f"Fragment found in {name}: {fragment}")
        else:
            print(f"No fragment found in {name} between {primer5} and {primer3}")


if __name__ == '__main__':
    main(
        fasta_file = "current_GREENGENES_gg16S_unaligned.fasta",
        primer_file="primers.fasta"
    )