# add header of sanger reads back to blast results

import pandas as pd

COL_NAMES = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


def read_tabular(path, header=None):
    if header:
        return pd.read_csv(str(path), sep="\t", names=header)
    else:
        return pd.read_csv(str(path), sep="\t")


def main():
    blast_table = read_tabular(snakemake.input["blast"], COL_NAMES)
    lookup_table = snakemake.params["samples"]

    output_table = str(snakemake.output)
    merge = pd.merge(blast_table, lookup_table, left_on="sseqid", right_on="hash")
    merge.to_csv(output_table, sep="\t", index=None)


if __name__ == "__main__":
    main()
