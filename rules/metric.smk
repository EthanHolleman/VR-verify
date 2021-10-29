
rule test_and_digest:
    # input blast alignments, filter by metrics in the config
    # file and simulate digestion with EcoRI and SacI
    conda:
        '../envs/Py.yml'
    input:
        blast='output/blast-results/with-headers/{insert_name}/{seq_file}/{insert_name}-{seq_file}.tsv',
        phred='output/metrics/phred/phred-quality-scores-all-reads.tsv'
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
       expand(
                'output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv',
                seq_file='all-runs-fasta-concat-md5', insert_name=inserts
                )
    output:
        'output/metrics/concat-metrics/concat-metrics.tsv'
    script:'../scripts/concat_with_headers.py'


rule make_global_summary_plots:
    conda:
        '../envs/R.yml'
    input:
        metrics_table='output/metrics/concat-metrics/concat-metrics.tsv',
        insert_list='output/VR-refs/all-insert-names.tsv'
    params:
        samples_table=config['samples_table']
    output:
        'output/metrics/global-summary/global-summary-plots.png'
    script:'../scripts/globalSummary.R'


rule extract_phred_scores_from_ab1_files:
    conda:
        '../envs/Py.yml'
    params:
        samples=SAMPLES
    output:
        'output/metrics/phred/phred-quality-scores-all-reads.tsv'
    script:'../scripts/extract_ab1_phred.py'


rule make_plots:
    conda:
        '../envs/R.yml'
    input:
        metrics_table='output/metrics/test-digest/{insert_name}/{seq_file}/{insert_name}-{seq_file}-filtered-blast-table-with-digests.tsv',
        phred_table='output/metrics/phred/phred-quality-scores-all-reads.tsv'
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

