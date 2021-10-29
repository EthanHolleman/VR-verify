from pathlib import Path
from Bio import SeqIO
import pandas as pd


RUNS = "runs/"
OUTPUT_TABLE = "samples/runs.tsv"

PRIMER_DISTS = {
    # Expected distance in nucleotides between end of primer binding site
    # and start of variable region insert for a
    # complete construct.
    "pFC9_t7_primer_1": 105,
    "Variable_region_insert_primer_1": 0,
    "Variable_region_insert_primer_2": 0,
    "VRI-1": 0,
}


def get_run_pairs(run_dir):
    files = [f for f in run_dir.iterdir() if f.stem != "all_text"]
    pairs = {}  # pairs of fasta and ab1 files
    pairs = {f.stem: [None, None] for f in files}
    for f in files:
        if f.suffix == ".seq":  # Fasta
            pairs[f.stem][0] = f
        elif f.suffix == ".ab1":
            pairs[f.stem][1] = f

    for each_stem, each_pair in pairs.items():
        assert each_pair != [None, None], each_stem
    return list(pairs.values())


def parse_run_dir_name(run_dir):
    # parse run dir names for run information
    date = Path(run_dir).name.split("_")[0]
    return {"date": date}


def parse_fasta_file(fasta_path):
    record = SeqIO.read(str(fasta_path), "fasta")
    read_name = record.description.split("_")[0].strip()
    primer = record.description.split("__")[0].replace(read_name, "")[1:].strip()
    return {
        "read_name": read_name,
        "primer": primer,
        "expected_insert_distance": PRIMER_DISTS[primer],
        "fasta_path": str(fasta_path),
        "hash": str(hash(str(fasta_path))),
        "read_length": len(record),
        "run_name": Path(fasta_path).parent.stem,
    }


def parse_abi_file(abi_path):
    record = SeqIO.read(str(abi_path), "abi")
    annos = record.annotations
    return {
        "sample_well": str(annos["sample_well"], "utf-8"),
        "dye": str(annos["dye"], "utf-8"),
        "run_start": annos["run_start"],
        "run_end": annos["run_finish"],
        "abi_path": str(abi_path),
    }


def file_pair_to_row(parent_dir, file_pair):

    fields = [
        parse_fasta_file(file_pair[0]),
        parse_abi_file(file_pair[1]),
        parse_run_dir_name(str(parent_dir)),
    ]
    row = {}
    for d in fields:
        row.update(d)
    return row


def parse_all_runs(runs_parent):
    rows = []
    for each_run in Path(runs_parent).iterdir():
        file_pairs = get_run_pairs(each_run)
        for each_file_pair in file_pairs:
            print(each_file_pair)
            rows.append(file_pair_to_row(each_run, each_file_pair))

    return pd.DataFrame(rows)


def main():

    parse_all_runs(RUNS).to_csv(OUTPUT_TABLE, sep="\t", index=None)


if __name__ == "__main__":
    main()
