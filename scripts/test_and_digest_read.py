from Bio.Restriction import *
from Bio.Restriction import EcoRI, SacI
from Bio.Seq import Seq
import pandas as pd
from Bio import SeqIO


def merge_mean_phred_score(table, phred_table):
    phred_means = phred_table.groupby("hash")["phred"].mean().reset_index()
    merge = pd.merge(table, phred_means, on="hash")
    print(merge)
    return merge


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


def make_decision(
    read_hash,
    mean_phred,
    primer,
    digest_length,
    start_diff,
    pident,
    align_length,
    expected_digest_length=223,
    max_digest_diff=10,
    min_pident=95,
    min_mean_phred=20,
    min_align_length=200,
):
    total_points = 0
    # Consider quality of read
    if mean_phred < min_mean_phred:
        total_points -= 1
    if primer == "pFC9_t7_primer_1":
        if abs(digest_length - expected_digest_length) > max_digest_diff:
            total_points -= 1
    if pident < min_pident:
        total_points -= 2
    if align_length < min_align_length:
        total_points -= 2

    status = ""
    if total_points == 0:
        status = "High confidence"
    elif total_points == -1:
        status = "Intermediate confidence"
    else:
        status = "Low confidence"

    return {"status": status, "points": total_points, "hash": read_hash}


def make_decision_all_reads(output_table):
    table_dict = output_table.to_dict(orient="records")
    verdict = []
    for er in table_dict:
        verdict.append(
            make_decision(
                er["hash"],
                er["phred"],
                er["primer"],
                er["EcoRI_SacI_digestion_frag_length"],
                er["length"],
                er["pident"],
                er["length"],
            )
        )

    # handle case were we have no alignments and therefore verdict_df made from empty list
    if not verdict:
        verdict_df = pd.DataFrame(columns=["hash", "status"])
    else:
        verdict_df = pd.DataFrame(verdict)

    return pd.merge(output_table, verdict_df, on="hash")


def get_candidate_fastas(filtered_blast_table):
    # Return paths to fasta files of Sanger reads that survived
    # filtering.
    return list(filtered_blast_table["fasta_path"])


def main():

    blast_alignments_merge = snakemake.input["blast"]

    phred_scores_table = snakemake.input["phred"]
    phred_scores_df = pd.read_csv(phred_scores_table, sep="\t")

    config = snakemake.params["config"]
    filtered_alignments = filter_alignments(
        blast_alignments_merge,
        config["blast_alignment_filters"]["min_length"],
        config["blast_alignment_filters"]["min_pident"],
    )
    print("=" * 10)
    print(len(filtered_alignments))
    filered_alignments_phred = merge_mean_phred_score(
        filtered_alignments, phred_scores_df
    )
    filtered_aligns_RE_digest = add_EcoRI_SacI_digest_column(filered_alignments_phred)
    filtered_aligns_RE_digest_verdict = make_decision_all_reads(
        filtered_aligns_RE_digest
    )

    output_table_path = str(snakemake.output)
    filtered_aligns_RE_digest_verdict.to_csv(output_table_path, sep="\t", index=None)


if __name__ == "__main__":
    main()
