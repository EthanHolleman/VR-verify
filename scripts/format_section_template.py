import pandas as pd
from pathlib import Path


def read_template(template_path):
    with open(str(template_path)) as handle:
        return handle.read()


def write_formatted_template(formated_template_string, output_path):
    with open(str(output_path), "w") as handle:
        handle.write(formated_template_string)
    return output_path


def format_template(template_string, format_values_dict):
    return template_string.format(**format_values_dict)


def summary_plot(summary_plot_path, output_path):

    p = str(Path("../" * 6).joinpath(summary_plot_path))
    return {"summary_plot_path": p}


def make_title(insert_name):
    return {"section_title": insert_name.strip(), "insert_name": insert_name.strip()}


def number_reads(blast_table):
    return {"number_of_reads": len(blast_table)}


def make_summary_table(blast_table):
    columns = blast_table[
        ["read_name", "run_name", "date", "primer", "read_length", "sample_well", "dye"]
    ]
    return {"summary_table": columns.to_markdown(index=False)}


def make_confidence_table(blast_table):
    columns = blast_table[["read_name", "date", "status"]]
    return {"confidence_table": columns.to_markdown(index=False)}


def main():

    output_path = str(snakemake.output)

    blast_table = pd.read_csv(snakemake.input["blast_table"], sep="\t")
    summary_plot_png = snakemake.input["summary_plot"]
    config = snakemake.params["config"]
    template_path = snakemake.params["template_path"]
    insert_name = snakemake.params["insert_name"]

    template_string = read_template(template_path)

    format_dict = {}
    format_dict.update(make_summary_table(blast_table))
    format_dict.update(make_title(insert_name))
    format_dict.update(number_reads(blast_table))
    format_dict.update(summary_plot(summary_plot_png, output_path))
    format_dict.update(make_confidence_table(blast_table))

    formated_template_string = format_template(template_string, format_dict)
    print(template_path.format(**format_dict))

    write_formatted_template(formated_template_string, output_path)


if __name__ == "__main__":
    main()
