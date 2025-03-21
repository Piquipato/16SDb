from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re, os, shutil
import json
import joblib as jb
from itertools import product
import tempfile
import pickle
import hashlib
import argparse
from datetime import datetime

"""
TO-DOs:
x Escribir una función que convierta una secuencia de nucleótidos con 
  ambigüedades a una expresión regular.
o Iterar sobre todos los primers presentes en el archivo primers.fasta.
o Iterar sobre todas las secuencias presentes en el archivo.
o Paralelizar la iteración sobre las secuencias del archivo.
o Guardar los fragmentos en un archivo fasta nuevo.
o Crear un pd.DataFrame con (bicho, región 16S) == fragmento.
o Guardar el DataFrame en un archivo csv/xlsx.
"""

def convert_ambiguous_to_regex(primer):
    REGEXER = json.load(open("pattern.json", "r"))
    regex = ""
    for i in range(len(primer)):
        regex += REGEXER[primer[i]]
    return regex



def find_fragment(sequence, primer5, primer3):
    primer5_regex = re.compile(convert_ambiguous_to_regex(primer5))
    primer3_regex = re.compile(convert_ambiguous_to_regex(primer3.reverse_complement()))

    match5 = primer5_regex.search(str(sequence))
    match3 = primer3_regex.search(str(sequence), match5.end() if match5 else 0)

    if match5 and match3:
        return Seq(sequence[match5.end():match3.start()])
    return None

def hash_sequence(sequence, size=10):
    h = hashlib.blake2b(digest_size=size)
    h.update(pickle.dumps(sequence))
    return h.hexdigest()

def worker(primer_regions, fasta_record, region, tempdir):
    primer5 = primer_regions[region]["F"]
    primer3 = primer_regions[region]["R"]
    seqhash = hash_sequence(fasta_record.seq)
    name = f"fragment_{fasta_record.id}_{fasta_record.name}_{seqhash}_{region}"
    description = "Region: {region} | ID: {id} | Name: {name} | Description: {description}".format(
        region=region,
        id=fasta_record.id,
        name=fasta_record.name,
        description=fasta_record.description
    )
    fragment = find_fragment(fasta_record.seq, primer5, primer3)
    if not isinstance(fragment, Seq) and fragment:
        fragment = Seq(fragment)
    fragment_file = None
    if fragment:
        fragment_record = SeqRecord(fragment, id=seqhash, name=name, description=description)
        fragment_file = os.path.abspath(f"{tempdir}/{name}.fasta")
        SeqIO.write(fragment_record, fragment_file, "fasta")


def main(fasta_file, primer_file, outfile="output.fasta", n_jobs=1, debug=False, logger=None):
    # Load fasta file
    if debug:
        print("Loading fasta files...")
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    primer_sequences = list(SeqIO.parse(open(primer_file), 'fasta'))
    
    # Organize primers by region
    primer_regions = {}
    for primer in primer_sequences:
        region, orientation = primer.name.split("_")
        if not region in primer_regions.keys():
            primer_regions[region] = {}
        primer_regions[region][orientation] = primer.seq
    
    tempdir = os.path.join(os.path.split(os.path.abspath(fasta_file))[0], "work/tmp")
    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
    os.mkdir(tempdir)

    # Find fragments for each region
    if debug:
        print("Starting Parallel Worker...")
    if n_jobs != 1:
        jb.Parallel(n_jobs=n_jobs, verbose=1)(
            jb.delayed(worker)(
                primer_regions, fasta_sequence, region, tempdir
            ) for fasta_sequence, region in product(fasta_sequences, list(primer_regions.keys()))
        )
    else:
        for fasta_sequence, region in product(fasta_sequences, list(primer_regions.keys())):
            worker(primer_regions, fasta_sequence, region, tempdir)

    if debug:
        print("Finished calculating fragments.")
        print("Merging fragments into a single file...")
    outfile = os.path.abspath(outfile)
    for fragment in os.listdir(tempdir):
        path = os.path.join(tempdir, fragment)
        temp_records = SeqIO.parse(open(path, "r"), 'fasta')
        for record in temp_records:
            SeqIO.write(record, open(outfile, "a"), "fasta")
        os.remove(path)
    shutil.rmtree(tempdir)

    if debug:
        print("Database created successfully.")
   

def cli():
    
    args = argparse.ArgumentParser()
    args.add_argument("--output", type=str, required=True, default="output.fasta")
    args.add_argument("--n_jobs", type=int, required=False, default=1)
    args.add_argument("--debug", action="store_true")

    args = args.parse_args()
    if args.debug:
        print("Creating database...")
        print(f"Output file: {args.output}")
        print(f"Number of jobs: {args.n_jobs}")
    main(
        fasta_file = "current_GREENGENES_gg16S_unaligned.fasta",
        primer_file="primers.fasta",
        outfile=args.output,
        n_jobs=args.n_jobs,
        debug=args.debug
    )



if __name__ == '__main__':
    cli()