import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

"""
This script performs a gene centric assembly using the gc-assembler of MEGAN of a certain SEED attribute e.g. Butyril-CoA
Dehydrogenae. It then takes the resulting contigs and uses NCBI-BLAST to find out the bacterial strains correlated to
the SEED attribute. The output is a list of organisms with the desired SEED attribute
"""

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('accession', type=str,
                    help='EC accession number of the desired SEED attribute for which the gc-assembly should be created')

args = parser.parse_args()

print(args.accession)


def gc_assembly(accession):
    
    return
