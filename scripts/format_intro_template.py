from datetime import datetime
from pathlib import Path
import pandas as pd
from format_section_template import (
    format_template,
    read_template,
    write_formatted_template,
)


def blast_params(config):
    blast_params = config["blast_alignment_filters"]
    return {
        "min_align_length": blast_params["min_length"],
        "min_pident": blast_params["min_pident"],
    }


def current_date():
    return {"current_date": str(datetime.today())}


def summary_plot(summary_plot_path):
    p = Path("../" * 3).joinpath(summary_plot_path)
    return {"summary_plot_png": p}


def confidence_summary_table(metrics_table, all_inserts_list):
    inserts = list(metrics_table["qseqid"])
    best_results = []
    for each_insert in all_inserts_list:
        reads = metrics_table.loc[metrics_table["qseqid"] == each_insert]
        if len(reads) > 0:
            reads = reads.sort_values(by="points", ascending=False)
            best_status = reads.iloc[0]
            best_results.append(
                {"Insert": each_insert, "Confidence": best_status["status"]}
            )
        else:
            best_results.append(
                {"Insert": each_insert, "Confidence": "No mapped reads"}
            )

    best_results = pd.DataFrame(best_results).sort_values(by="Confidence")
    # use dataframe just created which shows best status for each insert
    # to get number of inserts per status table "completeness"

    confidence_levels = [
        "High confidence",
        "Intermediate confidence",
        "Low confidence",
        "No mapped reads",
    ]
    completeness = [
        {
            "Number reads": len(best_results.loc[best_results["Confidence"] == cl]),
            "Confidence": cl,
        }
        for cl in confidence_levels
    ]

    return {
        "confidence_table": best_results.to_markdown(index=False),
        "completeness_table": pd.DataFrame(completeness).to_markdown(index=False),
        "percent_complete": round(
            len(best_results.loc[best_results["Confidence"] == confidence_levels[0]])
            / len(all_inserts_list),
            3,
        )
        * 100,
    }


def main():

    config = snakemake.params["config"]
    template = snakemake.params["intro_template"]
    summary_plot_path = snakemake.input["summary_plot"]
    metrics = pd.read_csv(snakemake.input["metrics"], sep="\t")
    inserts = snakemake.params["all_inserts"]

    output_path = str(snakemake.output)
    template_string = read_template(template)

    format_dict = {}
    format_dict.update({"title": config["title"]})
    format_dict.update(blast_params(config))
    format_dict.update(current_date())
    format_dict.update(summary_plot(summary_plot_path))
    format_dict.update(confidence_summary_table(metrics, inserts))

    formatted_template_string = format_template(template_string, format_dict)
    write_formatted_template(formatted_template_string, output_path)


if __name__ == "__main__":
    main()
