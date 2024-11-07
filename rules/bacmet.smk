rule bacmet:
    input:
        gene_proteins = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_proteins.faa",
    output:
        bacmet_result = f"{PROJECTNAME}/{{sample}}/annotation/bacmet/{{sample}}_BRG_MRG_bacmet.txt"
    params:
        database_path = config["bacmet"]["database_path"],
        bacmet_scan = config["bacmet"]["bacmet_scan"]
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_bacmet.log"
    conda:
        "../envs/blast.yaml"
    threads:
        config["bacmet"]["threads"]
    shell:
        """
        chmod +x {params.bacmet_scan}
        {params.bacmet_scan} -cpu {threads} -i {input.gene_proteins} -protein -blast -table -counts -o {output.bacmet_result} -d {params.database_path} 2>{log}
        """