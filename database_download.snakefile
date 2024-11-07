configfile: "config/config_database.json"

if not os.path.exists("database"):
    os.makedirs("database")
if not os.path.exists("database/logs"):
    os.makedirs("database/logs")


rule all:
    input:
        "metaphlan_db_downloaded",
        "bacmet_db_downloaded"

onsuccess:
    shell("rm *_db_downloaded")


rule metaphlan_database:
    output:
        "metaphlan_db_downloaded"
    log:
        "database/logs/metaphlan_database.log"
    conda:
        "envs/metaphlan.yaml"
    params:
        database_path = config['metaphlan_database']['database_path']
    shell:
        """
        rm -rf {params.database_path}
        mkdir database/metaphlan_database
        metaphlan --install --bowtie2db {params.database_path} 2>{log}
        touch metaphlan_db_downloaded
        """

rule bacmet_database:
    output:
        "bacmet_db_downloaded"
    params:
        database_path = config['bacmet_database']['database_path'],
        link_fasta = config['bacmet_database']['path_fasta'],
        link_annotation = config['bacmet_database']['path_annotation'],
        link_bacmet_scan = config['bacmet_database']['path_bacmet_scan']
    shell:
        """
        mkdir {params.database_path}
        wget -P {params.database_path} {params.link_fasta}
        wget -P {params.database_path} {params.link_annotation}
        wget -P {params.database_path} {params.link_bacmet_scan}
        touch bacmet_db_downloaded
        """
    
        