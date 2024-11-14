rule metaquast:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz",
        clean_metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_filtered_contig.fa"
    output:
        report = f"{PROJECTNAME}/{{sample}}/metagenome/metaquast/report.html"
    params:
        output_folder = f"{PROJECTNAME}/{{sample}}/metagenome/metaquast/"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_metaquast.log"
    conda:
        "../envs/quast.yaml"
    singularity:
        "docker://quay.io/biocontainers/quast:5.2.0--py38pl5321h40d3509_4"
    threads:
        config['metaquast']['threads']
    shell:
        """
        metaquast -t {threads} \
        -1 {input.clean_host_R1} \
        -2 {input.clean_host_R2} \
        -o {params.output_folder} \
        {input.clean_metagenome} \
        2>{log}
        """