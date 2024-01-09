rule regroup_file_metaphlan:
    input:
        metaphlan_result = expand(f"{PROJECTNAME}/{{sample}}/metaphlan/{{sample}}_metaphlan_profile.txt", sample=SAMPLES)
    output:
        txtfile_path = f"{PROJECTNAME}/sample_analysis/bacterial_community/file_path_metaphlan.txt"
    shell:
        """
        realpath {input.metaphlan_result} >> {output.txtfile_path}
        """

rule bacterial_population:
    input:
        txtfile_path = f"{PROJECTNAME}/sample_analysis/bacterial_community/file_path_metaphlan.txt"
    output:
        result_html = f"{PROJECTNAME}/sample_analysis/bacterial_community/result_ecological_analysis.html"
    params:
        threshold_best_abundance = config['r-ecological']['threshold_best_abundance'],
        sample_directory = f"{PROJECTNAME}/sample_analysis/bacterial_community/",
        project_directory = f"{PROJECTNAME}"
    log:
        f"{PROJECTNAME}/logs/ecological_analysis.log"
    conda:
        "../envs/r-ecological.yaml"
    shell:
        """
        Rscript -e "rmarkdown::render('scripts/ecological_analyses.Rmd', output_dir={params.sample_directory}, output_file='result_ecological_analysis.html', params=list(directory = {params.sample_directory}, threshold_best_abundance = {params.threshold_best_abundance}))"
        """