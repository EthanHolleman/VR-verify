from Bio import SeqIO
import pandas as pd


def extract_phred(ab1_path, samples):
    record = SeqIO.read(ab1_path, "abi")
    record_hash = int(samples.loc[samples["abi_path"] == ab1_path]["hash"])
    return pd.DataFrame(
        [
            {"phred": phred, "position": i + 1, "hash": record_hash}
            for i, phred in enumerate(record.letter_annotations["phred_quality"])
        ]
    )


def main():

    samples = snakemake.params["samples"]
    abi_files = list(samples["abi_path"])
    phred_scores = [extract_phred(abi, samples) for abi in abi_files]
    phred_table = pd.concat(phred_scores)
    phred_table.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()
