rule metagenome_bowtie_build:
    input:
        metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome.contigs.fa"
    output:
        index_metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome.1.bt2"
    params:
        metagenome_basename = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_metagenome_bowtie_build.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4"
    threads:
        config['metagenome_bowtie_build']['threads']
    shell:
        """
        bowtie2-build --threads {threads} {input.metagenome} {params.metagenome_basename} 2>{log}
        """

rule alignment_metagenome:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz",
        metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome.contigs.fa"
    output:
        metagenome_coverage = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome.cov"
    params:
        metagenome_basename = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_metagenome_coverage.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:4a9a541c36b8cd94820e79a67718c44f3522b7c6-0"
    threads:
        config['metagenome_coverage']['threads']
    shell:
        """
        bowtie2 -p {threads} -x {params.metagenome_basename} -1 {input.clean_host_R1} -2 {input.clean_host_R2} | \
        samtools view -b -h -@ {threads} - | \
        samtools sort -@ {threads} - | \
        samtools coverage --reference {input.metagenome} -o {output.metagenome_coverage} - 2>>{log}
        """

rule filter_contig:
    input:
        metagenome_coverage = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome.cov"
    output:
        contig_list = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_contig_ID.list"
    params:
        mean_coverage = config['filter_contig']['mean_coverage'],
        contig_len = config['filter_contig']['contig_len']
    shell:
        """
        awk '$7>={params.mean_coverage} && $3>={params.contig_len}' {input.metagenome_coverage} | \
        awk '{{print $1}}' > {output.contig_list}
        """

rule clean_metagenome:
    input:
        contig_list = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_contig_ID.list",
        metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome.contigs.fa"
    output:
        clean_metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_filtered_contig.fa"
    conda:
        "../envs/bowtie2.yaml"
    container:
        "docker://quay.io/biocontainers/seqtk:1.4--he4a0461_2"
    shell:
        """
        seqtk subseq {input.metagenome} {input.contig_list} > {output.clean_metagenome}
        """