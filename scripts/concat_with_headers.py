import pandas as pd


def main():

    input_files = snakemake.input
    dfs = [pd.read_csv(f, sep="\t") for f in input_files]
    concat = pd.concat(dfs)
    concat.to_csv(str(snakemake.output), sep="\t", index=None)


if __name__ == "__main__":
    main()
