rule metaquast:
    input:
        clean_metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_filtered_contig.fa"
    output:
        report = f"{PROJECTNAME}/{{sample}}/metagenome/metaquast/report.html"
    params:
        output_folder = f"{PROJECTNAME}/{{sample}}/metagenome/metaquast/"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_metaquast.log"
    conda:
        "../envs/quast.yaml"
    threads:
        config['metaquast']['threads']
    shell:
        """
        metaquast -t {threads} -o {params.output_folder} {input.clean_metagenome} 2>{log}
        """