from Bio.Restriction import *
from Bio.Restriction import EcoRI, SacI
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO


def EcoRI_SacI_digest(fasta_path):
    record = SeqIO.read(fasta_path, "fasta")
    rb = RestrictionBatch([EcoRI, SacI])

    sites = rb.search(record.seq)  # find restriction sites

    digest_dict = {"fasta_path": fasta_path}
    if len(sites[EcoRI]) == 1 and len(sites[SacI]) == 1:
        digest_dict.update(
            {"EcoRI_SacI_digestion_frag_length": abs(sites[EcoRI][0] - sites[SacI][0])}
        )
    else:
        digest_dict.update(
            {"EcoRI_SacI_digestion_frag_length": -1}
        )  # failed to digest in expected manor
    return digest_dict


def add_EcoRI_SacI_digest_column(filtered_blast_table):
    fasta_paths = get_candidate_fastas(filtered_blast_table)
    digests = [EcoRI_SacI_digest(fp) for fp in fasta_paths]
    if len(digests) == 0:
        digests_table = pd.DataFrame(
            columns=["fasta_path", "EcoRI_SacI_digestion_frag_length"]
        )
    else:
        digests_table = pd.DataFrame(digests)

    return pd.merge(filtered_blast_table, digests_table, on="fasta_path")


def filter_alignments(blast_table, min_length, min_pident):
    # Filter all alignments of sanger sequences -> VR reference by
    # min length and percent identity.
    blast_df = pd.read_csv(blast_table, sep="\t")
    return blast_df.loc[
        (blast_df["pident"] >= min_pident) & (blast_df["length"] >= min_length)
    ]


def get_candidate_fastas(filtered_blast_table):
    # Return paths to fasta files of Sanger reads that survived
    # filtering.
    return list(filtered_blast_table["fasta_path"])


def main():

    blast_alignments_merge = snakemake.input["blast"]
    config = snakemake.params["config"]
    filtered_alignments = filter_alignments(
        blast_alignments_merge,
        config["blast_alignment_filters"]["min_length"],
        config["blast_alignment_filters"]["min_pident"],
    )
    filtered_aligns_RE_digest = add_EcoRI_SacI_digest_column(filtered_alignments)
    output_table = str(snakemake.output)
    filtered_aligns_RE_digest.to_csv(output_table, sep="\t", index=None)


if __name__ == "__main__":
    main()
