rule alignement_host_genome:
    input:
        index_host_genome = f"{PROJECTNAME}/host_genome/host_genome.1.bt2",
        clean_R1 = f"{PROJECTNAME}/{{sample}}/reads/fastp/{{sample}}_clean_R1.fastq.gz",
        clean_R2 = f"{PROJECTNAME}/{{sample}}/reads/fastp/{{sample}}_clean_R2.fastq.gz"
    output:
        alignement_host_bam = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_reads.bam"
    params:
        host_basename = f"{PROJECTNAME}/host_genome/host_genome"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_alignement_host_genome.log"
    conda:
        "../envs/bowtie2.yaml"
    container:
        "docker://quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:4a9a541c36b8cd94820e79a67718c44f3522b7c6-0"
    threads:
        config['alignement_host_genome']['threads']
    shell:
        """
        bowtie2 -p {threads} \
        -x {params.host_basename} \
        -1 {input.clean_R1} \
        -2 {input.clean_R2} | \
        samtools view -b -@ {threads} - > {output.alignement_host_bam} \
        2>{log}
        """