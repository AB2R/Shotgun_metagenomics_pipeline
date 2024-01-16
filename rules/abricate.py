rule abricate:
    input:
        gene_nucleotides = f"{PROJECTNAME}/{{sample}}/annotation/prodigal/{{sample}}_genes_nucleotides.fna"
    output:
        card_ARG = f"{PROJECTNAME}/{{sample}}/annotation/abricate/{{sample}}_ARG.card.tab",
        ncbi_ARG = f"{PROJECTNAME}/{{sample}}/annotation/abricate/{{sample}}_ARG.ncbi.tab",
        resfinder_ARG = f"{PROJECTNAME}/{{sample}}/annotation/abricate/{{sample}}_ARG.resfinder.tab"
    params:
        min_identity = config['abricate']['minid'],
        min_coverage = config['abricate']['mincov']
    log:
        f"{PROJECTNAME}/logs/{{sample}}/{{sample}}_abricate.log"
    conda:
        "../envs/abricate.yaml"
    threads:
        config['abricate']['threads']
    shell:
        """
        abricate --threads {threads} --minid {params.min_identity} --mincov {params.min_coverage} --nopath --db card {input.gene_nucleotides} > {output.card_ARG}
        abricate --threads {threads} --minid {params.min_identity} --mincov {params.min_coverage} --nopath --db ncbi {input.gene_nucleotides} > {output.ncbi_ARG}
        abricate --threads {threads} --minid {params.min_identity} --mincov {params.min_coverage} --nopath --db resfinder {input.gene_nucleotides} > {output.resfinder_ARG}
        """