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
    threads:
        config['genes_bowtie_build']['threads']
    shell:
        """
        bowtie2-build --threads {threads} {input.genes} {params.genes_basename} 2>{log}
        """

rule reads_sort:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz"
    output:
        reads_sort_bam = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam"
    params:
        genes_basename = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_reads_abundance.log"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config['abundance']['threads']
    shell:
        """
        bowtie2 -p {threads} -x {params.metagenome_basename} -1 {input.clean_host_R1} -2 {input.clean_host_R2} | \
        samtools view -b -h -@ {threads} - | \
        samtools sort -@ {threads} -o {output.reads_sort_bam} - 2>{log}
        """

rule reads_abundance:
    input:
        reads_sort_bam = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam"
    output:
        reads_abundance = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_reads_abundance.tab"
        bam_index = f"{PROJECTNAME}/{{sample}}/annotation/abundance/{{sample}}_sort_reads.bam.bai"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_reads_abundance.log"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config['abundance']['threads']
    shell:
        """
        samtools index -@ {threads} -o {output.bam_index} {input.reads_sort_bam} 2>>{log}
        samtools idxstats -@ {threads} {input.reads_sort_bam} > {output.reads_abundance} 2>>{log}
        """