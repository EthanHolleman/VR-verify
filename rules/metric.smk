
rule test_and_digest:
    # input blast alignments, filter by metrics in the config
    # file and simulate digestion with EcoRI and SacI
    conda:
        '../envs/Py.yml'
    input:
        blast='output/blast-results/with-headers/{insert_name}/{seq_file}/{insert_name}-{seq_file}.tsv'
    output:
        'output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv'
    params:
        config=config
    script:
        '../scripts/test_and_digest_read.py'


rule merge_blast_digest_files:
    conda:
        '../envs/Py.yml'
    input:
        dynamic(expand(
                'output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv',
                seq_file='all-runs-fasta-concat-md5', allow_missing=True
                )
        )
    output:
        'output/metrics/concat-metrics/concat-metrics.tsv'
    script:'../scripts/concat_with_headers.py'


rule make_plots:
    conda:
        '../envs/R.yml'
    input:
        metrics_table='output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv'
    output:
        png='output/metrics/plots/{insert_name}/{seq_file}/{insert_name}-{seq_file}-plots.png',
        pdf='output/metrics/plots/{insert_name}/{seq_file}/{insert_name}-{seq_file}-plots.pdf'
    script:'../scripts/metricPlot.R'


rule make_all_metric_plots:
    input:
         dynamic(expand(
                'output/metrics/plots/{insert_name}/{seq_file}/{insert_name}-{seq_file}-plots.png',
                seq_file='all-runs-fasta-concat-md5', allow_missing=True
                )
        )
    output:
        'output/metrics/plots/.done.txt'
    shell:'''
    touch {output}
    '''

