#!/usr/bin/python3

import os
import sys
import argparse


def get_parser():
    """
    Initialise argument of Python script.
    
    Return:
        parser: Argument object who stock all information about argument.
    """
    parser = argparse.ArgumentParser(description="Download NCBI RefSeq genome.")

    parser.add_argument('-i', action='store', dest='id', required=True, type=str, help="NCBI RefSeq")

    return parser

def main():

    parser=get_parser()

    #Check if the arguments are present
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #Get RefSeq NCBI ID
    Args=parser.parse_args()
    refseq = Args.id

    os.chdir('host_genome')

    os.system(f'datasets download genome accession {refseq} --include genome')
    os.system(f'unzip ncbi_dataset.zip')
    os.system(f'find ncbi_dataset -name "*_genomic.fna" -exec mv -t . {{}} +')
    os.system(f'mv *genomic.fna {refseq}_genome.fna')
    os.system(f'rm -rf ncbi_dataset README.md')

if __name__ == "__main__":
    main()