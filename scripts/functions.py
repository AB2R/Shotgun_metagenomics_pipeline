import os


def mkdirectory(project_name, samples_list):
    """
    Make directory used for the pipeline.

    Args:
        project_name (str): The name of the project.
        samples_list (list): List of samples input in the pipeline.
    """
    if not os.path.exists(project_name):
        os.makedirs(project_name)

    list_project_directory = ['logs', 'host_genome', 'sample_analysis']
    for directory in list_project_directory:
        if directory == 'sample_analysis':
            if not os.path.exists(f"{project_name}/sample_analysis/bacterial_community"):
                os.makedirs(f"{project_name}/sample_analysis/bacterial_community")
            if not os.path.exists(f"{project_name}/sample_analysis/AMR"):
                os.makedirs(f"{project_name}/sample_analysis/AMR")

        elif directory == 'host_genome':
            if not os.path.exists(f"{project_name}/host_genome"):
                os.makedirs(f"{project_name}/host_genome")

    list_sample_directory = ['reads', 'metaphlan', 'metagenome', 'annotation']
    list_reads_directory = ['fastp', 'cleaning_host']
    list_metagenome_directory = ['megahit', 'clean_metagenome', 'metaquast']
    list_annotation_directory = ['prodigal', 'abricate', 'bacmet', 'abundance']
    for sample in samples_list:
        if not os.path.exists(f"{project_name}/logs/{sample['sampleID']}"):
            os.makedirs(f"{project_name}/logs/{sample['sampleID']}")

        for directory in list_sample_directory:

            if directory == 'reads':
                for d in list_reads_directory:
                    if not os.path.exists(f"{project_name}/{sample['sampleID']}/reads/{d}"):
                        os.makedirs(f"{project_name}/{sample['sampleID']}/reads/{d}")
            
            elif directory == 'metaphlan':
                if not os.path.exists(f"{project_name}/{sample['sampleID']}/metaphlan"):
                        os.makedirs(f"{project_name}/{sample['sampleID']}/metaphlan")
            
            elif directory == 'metagenome':
                for d in list_metagenome_directory:
                    if not os.path.exists(f"{project_name}/{sample['sampleID']}/metagenome/{d}"):
                        os.makedirs(f"{project_name}/{sample['sampleID']}/metagenome/{d}")

            elif directory == 'annotation':
                for d in list_annotation_directory:
                    if not os.path.exists(f"{project_name}/{sample['sampleID']}/annotation/{d}"):
                        os.makedirs(f"{project_name}/{sample['sampleID']}/annotation/{d}")



def get_number_sample(samples_list):
    """
    Return the number of samples input in the pipeline.

    Args:
        samples_list (list): List of samples input in the pipeline.

    Return:
        n_sample (int): List of n element in samples_list.
    """
    return len(samples_list)

def get_list_sample(samples_list):
    """
    Return a list of samples id.

    Args:
        samples_list (list): List of samples input in the pipeline.
    Return:
        list_id_sample (list): List of samples id.
    """
    id_sample = []
    for sample in samples_list:
        id_sample.append(sample['sampleID'])
    return id_sample
    
def get_all_output_files(config):
    """
    Return a list of files for rule all of the snakemake pipeline.

    Args:
        config (json file): config file of the analysis.
    Return:
        list_files (list): list of files output by the pipeline.
    """
    
    list_files = []

    for sample in config['sample']:
        #reads_filtered
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/reads/cleaning_host/{sample['sampleID']}_clean_host_reads_R1.fastq.gz")
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/reads/cleaning_host/{sample['sampleID']}_clean_host_reads_R1.fastq.gz")
        
        #metaphlan
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/metaphlan/{sample['sampleID']}_metaphlan_profile.txt")
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/metaphlan/{sample['sampleID']}_metaphlan.bowtie2.bz2")

        #metagenome
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/metagenome/clean_metagenome/{sample['sampleID']}_metagenome_filtered_contig.fa")

        #metaquast
        if config['metaquast']['process_analysis'] == "True":
            list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/metagenome/metaquast/report.html")

        #gene annotation
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/prodigal/{sample['sampleID']}_genes_proteins.faa")
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/prodigal/{sample['sampleID']}_genes_nucleotides.fna")

        #reads abundance
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/abundance/{sample['sampleID']}_reads_abundance.tab")

        #abricate
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/abricate/{sample['sampleID']}_ARG.card.tab")
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/abricate/{sample['sampleID']}_ARG.ncbi.tab")
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/abricate/{sample['sampleID']}_ARG.resfinder.tab")

        #bacmet
        list_files.append(f"{PROJECTNAME}/{sample['sampleID']}/annotation/bacmet/{sample['sampleID']}_BRG_MRG_bacmet.txt")

    if config['r-ecological']['process_analysis'] == "True":
        list_files.append(f"{PROJECTNAME}/sample_analysis/bacterial_community/result_ecological_analysis.html")
    
    return list_files