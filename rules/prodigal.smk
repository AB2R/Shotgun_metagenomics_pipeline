rule prodigal:
    input:
        clean_metagenome = f"{PROJECTNAME}/{{sample}}/metagenome/clean_metagenome/{{sample}}_metagenome_filtered_contig.fa"
    output:
        gene_coords = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_coords.gbk",
        gene_proteins = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_proteins.faa",
        gene_nucleotides = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_nucleotides.fna"
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_prodigal.log"
    conda:
        "../envs/prodigal.yaml"
    container:
        "docker://biocontainers/prodigal:v1-2.6.3-4-deb_cv1"
    shell:
        """
        prodigal -i {input.clean_metagenome} \
        -o {output.gene_coords} \
        -a {output.gene_proteins} \
        -d {output.gene_nucleotides} \
        -p meta
        """