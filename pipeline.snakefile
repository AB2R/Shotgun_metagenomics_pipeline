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
include: "rules/fastp.py"
include: "rules/download_genome_NCBI.py"
include: "rules/host_bowtie_build.py"
include: "rules/alignement_host_genome.py"
include: "rules/remove_host_reads.py"
include: "rules/metaphlan.py"
include: "rules/megahit.py"
include: "rules/filter_metagenome.py"
include: "rules/metaquast.py"
include: "rules/prodigal.py"
include: "rules/abundance.py"
include: "rules/abricate.py"
include: "rules/bacmet.py"

list_files_output = get_all_output_files(config)

rule all:
    input:
        expand(f"{PROJECTNAME}/{{sample}}/reads/fastp/{{sample}}_fastp.{{extension}}", sample=SAMPLES, extension=["json", "html"]),
        list_files_output
