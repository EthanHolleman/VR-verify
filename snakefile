import pandas as pd

configfile: 'config.yml'


SAMPLES = pd.read_csv('samples/runs.tsv', sep='\t')

BLAST_SUFFI = [
    'nos'
]

wildcard_constraints:
   seqfile = '\w+'

include: 'rules/blast.smk'
include: 'rules/metric.smk'


rule all:
    input:
        'output/metrics/plots/.done.txt'
