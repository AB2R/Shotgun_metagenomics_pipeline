#!/bin/bash

# Help function
Help()
{
    echo "Download Refseq NCBI genome of host animal."
    echo
    echo "Usage: download_host_genome.sh in.path_host_genome in.refseq"
    echo
}

# Assign variable
Directory=$1
Refseq=$2

if [ -z $Directory ]
then
    echo "No directory"
    Help
    exit
fi

if [ -z $Refseq ]
then
    echo "No Refseq NCBI accession number."
    Help
    exit
fi

# Download genome
cd ${Directory}
datasets download genome accession ${Refseq} --include genome
unzip ncbi_dataset.zip
find ncbi_dataset -name *_genomic.fna -exec mv -t . {} +
mv *genomic.fna host_genome.fna
rm -rf ncbi_dataset README.md ncbi_dataset.zip
