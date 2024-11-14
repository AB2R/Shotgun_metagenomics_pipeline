def get_host_genome(config):
    if config['host']['host_path_file'] == "":
        return(f"{PROJECTNAME}/host_genome/host_genome.fna")
    else:
        return(config['host']['host_path_file'])

rule host_bowtie_build:
    input:
        host_genome = get_host_genome(config)
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