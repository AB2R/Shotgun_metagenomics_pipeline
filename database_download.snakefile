configfile: "config/config_database.json"

if not os.path.exists("database"):
    os.makedirs("database")
if not os.path.exists("database/logs"):
    os.makedirs("database/logs")


rule all:
    input:
        "metaphlan_db_downloaded"


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
    
        