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
from Bio.Blast import NCBIXML
import Bio.Blast.NCBIWWW

import matplotlib.pyplot as plt
import numpy as np

E_VALUE_THRESH = 0.04
GC_OUTPUT = "gc_output_temp"

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


"""
Runs NCBI BLAST for the output of gc-assembly vs NR database
"""


def blast(input):
    print("INPUT: " + str(input))
    for record in SeqIO.parse(input, "fasta"):
        result_handle = Bio.Blast.NCBIWWW.qblast("blastn", "nr", record.seq)
        blast_record = NCBIXML.read(result_handle)
        with open(input[:-6] + "_blast_results", "w") as writer:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        writer.writelines("****Alignment****")
                        writer.writelines("sequence:" + str(alignment.title))
                        writer.writelines("length:" + str(alignment.length))
                        writer.writelines("e value:" + str(hsp.expect))
                        writer.writelines(hsp.query[0:75] + "...")
                        writer.writelines(hsp.match[0:75] + "...")
                        writer.writelines(hsp.sbjct[0:75] + "...")
                        writer.write("\n")
                        writer.write("\n")
        break


# Run the gene centric assembly using the command line parameters
gc_assembly_file_name = GC_OUTPUT + "_ID_" + args.id + ".fasta"
gc_assembly(args.input, gc_assembly_file_name, args.id)

# Run blast with the gc assembly output vs nr
#blast(gc_assembly_file_name)
blast("reads-Acyl-CoA_dehydrogenase__short-chain_specific__EC_1.3.99.2.fasta")
