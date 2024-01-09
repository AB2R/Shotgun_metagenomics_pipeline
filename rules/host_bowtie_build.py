rule host_bowtie_build:
    input:
        host_genome = f"{PROJECTNAME}/host_genome/host_genome.fna"
    output:
        index_host_genome = f"{PROJECTNAME}/host_genome/host_genome.1.bt2"
    params:
        host_basename = f"{PROJECTNAME}/host_genome/host_genome"
    log:
        f"{PROJECTNAME}/logs/host_bowtie_build.log"
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config['host_bowtie_build']['threads']
    shell:
        """
        bowtie2-build --threads {threads} {input.host_genome} {params.host_basename} 2>{log}
        """