rule genes_bowtie_build:
    input:
        genes = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_nucleotides.fna"
    output:
        index_genes = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes.1.bt2"
    params:
        genes_basename = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_genes_bowtie_build.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4"
    threads:
        config['genes_bowtie_build']['threads']
    shell:
        """
        bowtie2-build --threads {threads} {input.genes} {params.genes_basename} 2>{log}
        """

rule reads_sort:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz",
        index_genes = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes.1.bt2"
    output:
        reads_sort_bam = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam"
    params:
        genes_basename = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_reads_abundance.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:4a9a541c36b8cd94820e79a67718c44f3522b7c6-0"
    threads:
        config['abundance']['threads']
    shell:
        """
        bowtie2 -p {threads} -x {params.genes_basename} -1 {input.clean_host_R1} -2 {input.clean_host_R2} | \
        samtools view -b -h -@ {threads} - | \
        samtools sort -@ {threads} -o {output.reads_sort_bam} - 2>{log}
        """

rule reads_abundance:
    input:
        reads_sort_bam = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam"
    output:
        reads_abundance = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_reads_abundance.tab"        
    params:
        bam_index = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam.bai",
        remove_file = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam*"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_reads_abundance.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    threads:
        config['abundance']['threads']
    shell:
        """
        samtools index -@ {threads} -o {params.bam_index} {input.reads_sort_bam} 2>>{log}
        samtools idxstats -@ {threads} {input.reads_sort_bam} > {output.reads_abundance} 2>>{log}
        rm {params.remove_file}
        """