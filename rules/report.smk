
rule make_section_markdown:
    conda:
        '../envs/Py.yml'
    input:
        blast_table='output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv',
        summary_plot='output/metrics/plots/{insert_name}/{seq_file}/{insert_name}-{seq_file}-plots.png',
    output:
        'output/report/markdown/sections/{insert_name}/{seq_file}/{insert_name}.md'
    params:
        config=config,
        template_path=config['markdown']['section_template'],
        insert_name=lambda wildcards: wildcards.insert_name
    script:'../scripts/format_section_template.py'


rule make_all_section_markdown:
    input:
       expand(
                'output/report/markdown/sections/{insert_name}/{seq_file}/{insert_name}.md',
                seq_file='all-runs-fasta-concat-md5', insert_name=inserts
            )
    output:
        'output/report/markdown/sections/.done.txt'
    shell:'''
    touch {output}
    '''


rule concat_all_section_files:
    conda:
        '../envs/Py.yml'
    input:
       expand(
                'output/report/markdown/sections/{insert_name}/{seq_file}/{insert_name}.md',
                seq_file='all-runs-fasta-concat-md5', insert_name=inserts
            )
    output:
        'output/report/markdown/all_sectons_concat.md'
    script:'../scripts/cat_markdown_ordered.py'


rule make_intro:
    conda:
        '../envs/Py.yml'
    input:
        summary_plot='output/metrics/global-summary/global-summary-plots.png',
        metrics='output/metrics/concat-metrics/concat-metrics.tsv'
    output:
        'output/report/markdown/intro.md'
    params:
        config=config,
        intro_template=config['markdown']['intro_template'],
        all_inserts=inserts
    script:'../scripts/format_intro_template.py'


rule make_complete_markdown_report:
    input:
        sections='output/report/markdown/all_sectons_concat.md',
        intro='output/report/markdown/intro.md'
    output:
        'output/report/markdown/complete_report.md'
    shell:'''
    cat {input.intro} {input.sections} > {output}
    '''


rule convert_complete_report_to_html:
    conda:
        '../envs/pandoc.yml'
    input:
        'output/report/markdown/complete_report.md'
    output:
        'output/report/html/complete_report.html'
    shell:'''
    pandoc -s {input} -o {output}
    '''

