"""
 *  Copyright (C) 2019 Timo Niklas Lucas
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

This script performs a gene centric assembly using the gc-assembler of MEGAN of a certain SEED attribute e.g. Butyril-CoA
Dehydrogenae. It then takes the resulting contigs and uses NCBI-BLAST to find out the bacterial strains correlated to
the SEED attribute. The output is a list of organisms with the desired SEED attribute
"""

import argparse
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline


import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('input', type=str,
                    help='Input file for the gene centric assembly (.daa file)')
parser.add_argument('id', type=str,
                    help='EC accession number of the desired SEED attribute for which the gc-assembly should be created')
args = parser.parse_args()


"""
Performs Gene Centric assembly for a specific SEED ID on a given .daa file
The function calls the gc-assembler tool from MEGAN
Returns name of the output file
"""


def gc_assembly(input_file: str, output_file: str, seed_id: str) -> object:
    command = "/Applications/MEGAN_/tools/gc-assembler --input " + input_file + \
              " -o " + output_file + \
              " -fun SEED -id " + seed_id
    os.system(command)
    print(command)
    return output_file

"""
Runs NCBI BLAST for the output of gc-assembly vs NR database
"""
def online_blast(input):
    result_handles = []
    for record in SeqIO.parse(input, "fasta"):
        print(record.seq)
        blastn_cline = NcbiblastnCommandline(query=input, db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")
        os.system(str(blastn_cline))

    return result_handles


# Run the gene centric assembly using the command line parameters
#gc_output = gc_assembly(args.input, args.input[:-5] + "_gc_assembly.txt", args.id)
results = online_blast("/Users/Timo/Desktop/Mappings/Mappings/output_test_gc_assembly.txt")
#print(results)

#handle = Bio.Blast.NCBIWWW.qblast("blastn", "nr", "ATGGATTTTCATCTGACAAATGAGCATCTGATGCTCCGCAAGATGTACCGCGAATTCGCCGAGACCGAAGTGAAGCCAATCGCCACTGAGATCGACGAAGAAGAGCGCTTCCCCATGGAGACCGTGGAGAAGATGGCCAAGCTGGGCATGATGGGCATTTACTTCCCNAAGGAGTATGGCGGCGCCGGCGGCGACGTCCTCTCNTATGCAATGTGCGTGGAGGAGCTGGCGAAGGTNTGCGGTACTACTGCCGTNGTCGTTTCCGCTCATACCTCCCTGTGCTGCGCTCCNATCTTTGAACATGGTACCGAGGAGCAGAAGAAAAAGTATCTGCCTGACCTGCTTTCCGGCAAGAAGATTGGCGCGTTCGGTCTGACCGAACCCGGNGCCGGAACGGATGCCTCNGGTCAGCAGACCATTGCCGTTCTGGAAGGTGACCATTACATCCTCAACGGTTCCAAGTGCTTCATTACCAACGGCAACGTGGCCGATACCTTCGTGGTATTCGCCATGACCGATAAGAAGCTGGGCAACCACGGCATTTCCGCTTTCATTGTTGAAAGGGATTTCCCCGGNTTCTCTCAGGGCAAGCATGAGCTGAAGATGGGCATCCGCGGCAGCTCTACCTGCGACCTGATTTTCGAGGACTGCATTGTACCGAAAGAGAACCTGCTGGGNAAAGAGGGCAAGGGCTTTAAGGTTGCCATGCAGACCCTGGACGGCGGNCGCATCGGCATTGCCTCCCAGGCNCTGGGNCTGGCCGAGGGTGCCATCGATGAAGCTGTCAAGTATGTCAAGGAGCGCGTNCAGTTCGGCCGCCCCCTGTCCAAGTTCCAGAACACTCAGTTCCAGCTGGCCGACATGCACTGCCGCACCCAGGCTGCGCAGTATCTGGTGTACGCCGCCGCCTGCAAGAAGCAGCTCCATGAGGATTACTCAATGGATGCGGCCATGGCCAAGCTGTTTGCCGCCGAGACTGCTTCTGATGTGACCCGCCGTGCCGTACAGCTGTTTGGCGGNTACGGTTACACCCGTGAGTATCCCGTGGAGCGCATGATGCGCGACGCNAAGATCACCGAGATCTACGAGGGTACTTCCGAAGTCCAGCGGATGGTTATCGCCAGCCACTTGGGN")
