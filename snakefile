import pandas as pd

configfile: 'config.yml'


SAMPLES = pd.read_csv('samples/runs.tsv', sep='\t')

BLAST_SUFFI = [
    'ndb'
]

inserts = list(pd.read_csv('inserts.tsv', sep='\t')['insert_name'])

print(inserts)

wildcard_constraints:
   seqfile = '\w+'

include: 'rules/blast.smk'
include: 'rules/metric.smk'
include: 'rules/report.smk'


rule all:
    input:
        'output/report/html/complete_report.html'
        
