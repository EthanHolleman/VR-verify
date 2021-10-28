# Max allowed local id in blast is 50
from os import rename
from Bio import SeqIO
import pandas as pd
import time


def rename_record(fasta_path, samples):
    record = SeqIO.read(fasta_path, "fasta")
    record_hash = samples.loc[samples["fasta_path"] == str(fasta_path)]["hash"].values[
        0
    ]
    descrip = record.description
    record.description = ""
    record.id = str(record_hash)
    return record


def rename_all_records(paths, samples):
    return [rename_record(path, samples) for path in paths]


def main():

    fasta = snakemake.input["fasta"]
    output_fasta = str(snakemake.output)
    samples = snakemake.params["samples"]

    hashed_records = rename_all_records(fasta, samples)
    SeqIO.write(hashed_records, output_fasta, "fasta")


if __name__ == "__main__":
    main()
