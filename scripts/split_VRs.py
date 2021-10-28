# split complete inserts fasta into individual records

from Bio import SeqIO
from pathlib import Path


def main():

    VR_refs = snakemake.input["VR_refs"]
    VR_records = SeqIO.parse(VR_refs, "fasta")
    output_dir = str(snakemake.params["output_dir"])
    for each_record in VR_records:
        name = each_record.description.split("insert_")[-1]
        each_record.description = ""
        each_record.id = name.strip()
        output_path = Path(output_dir).joinpath(name).with_suffix(".fa")
        SeqIO.write([each_record], output_path, "fasta")


if __name__ == "__main__":
    main()
