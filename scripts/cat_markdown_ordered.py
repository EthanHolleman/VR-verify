# Concat markdown files in order of VR insert

from pathlib import Path


def main():

    markdown_files = snakemake.input
    sorted_markdown = sorted(
        markdown_files, key=lambda p: int(Path(p).stem.split("VR-")[-1])
    )
    concat_markdown = ""
    for each_file in sorted_markdown:
        with open(each_file) as handle:
            concat_markdown += handle.read()

    with open(str(snakemake.output), "w") as handle:
        handle.write(concat_markdown)


if __name__ == "__main__":
    main()
