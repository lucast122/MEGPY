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
import ntpath
import os
import signal

import Bio
import Bio.Blast.NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML


class TimeoutException(Exception):  # Custom exception class
    pass


def timeout_handler(signum, frame):  # Custom signal handler
    raise TimeoutException

E_VALUE_THRESH = 0.04

GC_OUTPUT = "gc_temp"

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('input', type=str,
                    help='Input file for the gene centric assembly (.daa file)')
parser.add_argument('id', type=str,
                    help='EC accession number of the desired SEED attribute for which the gc-assembly should be created')
args = parser.parse_args()

MAPPING_FILE_PATH = args.input
MAPPING_FILE_NAME = ntpath.basename(MAPPING_FILE_PATH)[:-5]


"""
Performs Gene Centric assembly for a specific SEED ID on a given .daa file
The function calls the gc-assembler tool from MEGAN
Returns name of the output file
"""


def gc_assembly(input_file: str, output_file: str, seed_id: str):
    command = "/Applications/MEGAN_/tools/gc-assembler --input " + input_file + \
              " -o " + output_file + \
              " -fun SEED -id " + seed_id
    os.system(command)
    print(command)


"""
Runs NCBI BLAST for the output of gc-assembly vs NR database
"""


def blast(input, max_time):
    strains = []
    BLAST_OUTPUT_NAME = input[:-6] + str(MAPPING_FILE_NAME) + "_blast_results"
    count = 0
    print("RUNNING BLAST WITH INPUT: " + str(input))
    try:
        os.remove(BLAST_OUTPUT_NAME)
        print("Blast results already exist for file " + input + " deleting results before computing new.")
    except FileNotFoundError:
        print("No previous results found")

    for record in SeqIO.parse(input, "fasta"):
        # Only BLAST contigs that came from more than 20 reads and have more than 50bp
        if len(record.seq) < 50 or int(record.description.split(" ")[2].split("=")[1]) < 20:
            continue

        print("Running BLAST for record ", record.id, record.seq)
        # Change the behavior of SIGALRM
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(max_time)
        # This try/except loop ensures that
        #   you'll catch TimeoutException when it's sent.
        try:
            result_handle = Bio.Blast.NCBIWWW.qblast("blastn", "nr", record.seq, megablast=True)
            with open("blast_output_" + MAPPING_FILE_NAME + args.id + ".xml", "w") as out_handle:
                out_handle.write(result_handle.read())

            result_handle.close()
        except TimeoutException:
            print("BLAST run aborted, because it took longer than " + str(max_time) + " seconds.")
            with open(BLAST_OUTPUT_NAME, "a") as writer:
                writer.writelines("BLAST failed. Runtime exceeded maximum time")

            continue  # continue the for loop if function A takes more than defined number of seconds
        else:
            # Reset the alarm
            signal.alarm(0)

        # Now write the BLAST output to a file and create another file with a list of strains
        blast_record = NCBIXML.read(result_handle)
        with open(BLAST_OUTPUT_NAME, "a") as writer:
            for alignment in blast_record.alignments:
                strains.append(alignment.title)
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
    return strains



# Run the gene centric assembly using the command line parameters
gc_assembly_file_name = GC_OUTPUT + "_ID_" + args.id + ".fasta"
# gc_assembly(args.input, gc_assembly_file_name, args.id)

# Run blast with the gc assembly output vs nr
strains = blast(gc_assembly_file_name, 60)
strains = list(dict.fromkeys(strains))

# Now write all unique strains associated with wanted MEGAN ID to a file
STRAINS_FILE_NAME = "strains_" + str(MAPPING_FILE_NAME) + args.id
with open(STRAINS_FILE_NAME, "w") as f:
    for strain in strains:
        f.writelines(strain)

# create summary file to informs about present strains as well as strains associated with reverse-beta-oxidation





# blast("reads-Acyl-CoA_dehydrogenase__short-chain_specific__EC_1.3.99.2.fasta")
# blast("reads-Butyryl-CoA_dehydrogenase__EC_1.3.99.2.fasta")
