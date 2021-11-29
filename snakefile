import pandas as pd
import os

configfile: 'config.yml'




SAMPLES = pd.read_csv('samples/runs.tsv', sep='\t')

BLAST_SUFFI = [
    'nhr'
]

inserts = list(pd.read_csv('inserts.tsv', sep='\t')['insert_name'])



wildcard_constraints:
   seqfile = '\w+'

include: 'rules/blast.smk'
include: 'rules/metric.smk'
include: 'rules/report.smk'


rule all:
    input:
        'output/report/html/complete_report.html'
        
