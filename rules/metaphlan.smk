rule metaphlan:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz"
    output:
        metaphlan_result = f"{PROJECTNAME}/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.txt",
        metaphlan_bowtie2 = f"{PROJECTNAME}/{{sample}}/metaphlan/{{sample}}_metaphlan.bowtie2.bz2"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_metaphlan.log"
    conda:
        "../envs/metaphlan.yaml"
    singularity:
        "docker://qbioturin/metaphlan4:0.3.2"
    threads:
        config['metaphlan']['threads']
    shell:
        """
        metaphlan {input.clean_host_R1},{input.clean_host_R2} \
        --bowtie2out {output.metaphlan_bowtie2} \
        --nproc {threads} \
        --input_type fastq \
        -o {output.metaphlan_result} \
        2>{log}
        """