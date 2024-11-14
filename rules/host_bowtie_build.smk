rule host_bowtie_build:
    input:
        if config['host']['host_path_file'] == "":
            host_genome = f"{PROJECTNAME}/host_genome/host_genome.fna"
        else:
            host_genome = config['host']['host_path_file']
    output:
        index_host_genome = f"{PROJECTNAME}/host_genome/host_genome.1.bt2"
    params:
        host_basename = f"{PROJECTNAME}/host_genome/host_genome"
    log:
        f"{PROJECTNAME}/logs/host_bowtie_build.log"
    conda:
        "../envs/bowtie2.yaml"
    singularity:
        "docker://quay.io/biocontainers/bowtie2:2.5.4--h7071971_4"
    threads:
        config['host_bowtie_build']['threads']
    shell:
        """
        bowtie2-build --threads {threads} {input.host_genome} {params.host_basename} 2>{log}
        """