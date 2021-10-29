from pathlib import Path


rule split_VR_refs:
    conda:
        '../envs/Py.yml'
    input:
        VR_refs='resources/VR-refs/complete_inserts.fa'
    output:
        expand('output/VR-refs/split/{insert_name}.fa', insert_name=inserts)
    params:
        output_dir='output/VR-refs/split'
    script:'../scripts/split_VRs.py'


rule make_inserts_list:
    input:
        expand('output/VR-refs/split/{insert_name}.fa', insert_name=inserts)
    output:
        'output/VR-refs/all-insert-names.tsv'
    params:
        split_dir='output/VR-refs/split'
    shell:'''
    echo "insert_name" > {output}
    ls -1 {params.split_dir} | sed -e 's/\.fa$//' >> {output}
    '''



rule concat_run_fasta_files:
    input:
        list(SAMPLES['fasta_path'])
    output:
        'output/runs/concat/all-runs-fasta-concat.fa'
    shell:'''
    cat {input} > {output}
    '''

rule hash_headers_for_blast_because_its_annoying:
    conda:
        '../envs/Py.yml'
    input:
        fasta=list(SAMPLES['fasta_path'])
    output:
        'output/runs/concat-hash/all-runs-fasta-concat-md5.fa'
    params:
        samples=SAMPLES
    script:'../scripts/hash_read_headers.py'


rule make_blast_db:
    conda:
        '../envs/blast.yml'
    input:
        seqs='output/runs/concat-hash/{seq_file}.fa'
    output:
        expand(
            'output/blast-dbs/{seq_file}/{seq_file}.{blast_suffi}', 
            allow_missing=True, blast_suffi=BLAST_SUFFI
            )
    params:
        db_dir=lambda wildcards: f'output/blast-dbs/{wildcards.seq_file}',
        db_name=lambda wildcards: f'output/blast-dbs/{wildcards.seq_file}/{wildcards.seq_file}',
        seq_file_name = lambda wildcards: wildcards.seq_file
    shell:'''
    mkdir -p {params.db_dir}
    makeblastdb -in {input.seqs} -parse_seqids -title "{params.seq_file_name}" -dbtype nucl -out {params.db_name}
    '''


rule run_blast:
    conda:
        '../envs/blast.yml'
    input:
        db=expand(
            'output/blast-dbs/{seq_file}/{seq_file}.{blast_suffi}', 
            allow_missing=True, blast_suffi=BLAST_SUFFI
            ),
        insert_ref='output/VR-refs/split/{insert_name}.fa'
    params:
        db_path=lambda wildcards: f'output/blast-dbs/{wildcards.seq_file}/{wildcards.seq_file}',
        output_dir='output/blast-results'
    output:
        'output/blast-results/{insert_name}/{seq_file}/{insert_name}-{seq_file}.tsv'
    shell:'''
    mkdir -p {params.output_dir}
    cat {input.insert_ref} | blastn -db {params.db_path} -perc_identity 0 -outfmt 6 > {output}
    '''


rule merge_header_back_in:
    conda:
        '../envs/Py.yml'
    input:
        blast='output/blast-results/{insert_name}/{seq_file}/{insert_name}-{seq_file}.tsv'
    params:
        samples=SAMPLES
    output:
        'output/blast-results/with-headers/{insert_name}/{seq_file}/{insert_name}-{seq_file}.tsv'
    script:'../scripts/merge_Blast_results.py'


