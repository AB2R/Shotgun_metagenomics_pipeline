### Import configfiles
configfile: "config/config.json" #Config file of the project with all the samples
configfile: "config/config_tools.json" #Config file of tools options and thresholds.

include: "scripts/functions.py"
PROJECTNAME = config['project_name']

### Create folders
mkdirectory(PROJECTNAME, config['sample'])

### Wildcards
SAMPLES = get_list_sample(config['sample'])

### Rules

### Import rules
include: "rules/fastp.smk"
include: "rules/download_genome_NCBI.smk"
include: "rules/host_bowtie_build.smk"
include: "rules/alignement_host_genome.smk"
include: "rules/remove_host_reads.smk"
include: "rules/metaphlan.smk"
include: "rules/bacterial_population.smk"
include: "rules/megahit.smk"
include: "rules/filter_metagenome.smk"
include: "rules/metaquast.smk"
include: "rules/prodigal.smk"
include: "rules/abundance.smk"
include: "rules/abricate.smk"
include: "rules/bacmet.smk"

list_files_output = get_all_output_files(config)

rule all:
    input:
        expand(f"{PROJECTNAME}/{{sample}}/reads/fastp/{{sample}}_fastp.{{extension}}", sample=SAMPLES, extension=["json", "html"]),
        list_files_output
