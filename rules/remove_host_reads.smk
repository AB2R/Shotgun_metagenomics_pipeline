rule remove_host_reads:
    input:
        alignement_host_bam = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_reads.bam"
    output:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_alignement_host_genome.log"
    conda:
        "../envs/samtools.yaml"
    container:
        "docker://quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
    threads:
        config['remove_host_reads']['threads']
    shell:
        """
        samtools view -b -@ {threads} -f 12 -F 256 \
        {input.alignement_host_bam} | \
        samtools sort -n -m 5G -@ {threads} - | \
        samtools fastq -@ {threads} - \
        -1 {output.clean_host_R1} \
        -2 {output.clean_host_R2} \
        -n 2>{log}
        """