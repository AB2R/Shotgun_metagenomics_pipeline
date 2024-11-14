rule megahit:
    input:
        clean_host_R1 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R1.fastq.gz",
        clean_host_R2 = f"{PROJECTNAME}/{{sample}}/reads/cleaning_host/{{sample}}_clean_host_reads_R2.fastq.gz"
    output:
        metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/{{sample}}_metagenome.contigs.fa"
    params:
        kmin = config['megahit']['k-min'],
        kstep = config['megahit']['k-step'],
        min_contig_len = config['megahit']['min-contig-len'],
        options = config['megahit']['others_options'],
        prefix = f"{{sample}}_metagenome",
        output_folder = f"{PROJECTNAME}/{{sample}}/metagenome/megahit/"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_megahit.log"
    conda:
        "../envs/megahit.yaml"
    singularity:
        "docker://quay.io/biocontainers/megahit:1.2.9--h43eeafb_5"
    threads:
        config['megahit']['threads']
    shell:
        """
        megahit -t {threads} \
        -1 {input.clean_host_R1} \
        -2 {input.clean_host_R2} \
        --k-min {params.kmin} \
        --k-step {params.kstep} \
        --min-contig-len {params.min_contig_len} \
        {params.options} \
        -o {params.output_folder} \
        --out-prefix {params.prefix} \
        2>{log}
        """