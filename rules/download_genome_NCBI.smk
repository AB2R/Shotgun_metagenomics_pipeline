rule download_genome_host:
    output:
        genome_host = f"{PROJECTNAME}/host_genome/host_genome.fna"
    params:
        refseq_ID = config['host_ncbi_refseq'],
        host_directory = f"{PROJECTNAME}/host_genome/"
    log:
        f"{PROJECTNAME}/logs/download_genome_host.log"
    conda:
        "../envs/ncbi-download.yaml"
    shell:
        """
        chmod +x scripts/download_host_genome.sh
        scripts/download_host_genome.sh {params.host_directory} {params.refseq_ID} 2>{log}
        """
