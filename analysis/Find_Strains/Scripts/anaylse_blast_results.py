#  """
#   *  Copyright (C) 2019 Timo Niklas Lucas
#   *
#   *  This program is free software: you can redistribute it and/or modify
#   *  it under the terms of the GNU General Public License as published by
#   *  the Free Software Foundation, either version 3 of the License, or
#   *  (at your option) any later version.
#   *
#   *  This program is distributed in the hope that it will be useful,
#   *  but WITHOUT ANY WARRANTY; without even the implied warranty of
#   *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   *  GNU General Public License for more details.
#   *
#   *  You should have received a copy of the GNU General Public License
#   *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  """

# THIS SCRIPT IS USED TO EXTRACT THE PROTEINS AFTER BLASTING A GC-ASSEMBLY AGAINST A DATABASE USING BLAST
# THE SCRIPTS TAKES BLASTX OUTPUT IN TABULAR FORMAT AND PRINTS INTERESTING PROTEINS THAT THE CONTIGS GOT MAPPED TO
import argparse
import os

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('--input', type=str,
                    help='Input file (blast output file)')
args = parser.parse_args()

INPUT_FILE = args.input
CONTIG_DICT = {}
contig_hit = ""
contig_id = ""
first_name_word_index = 0
last_name_word_index = 0
protein_name = ""
identity = ""


def gc_assembly(input_file: str, output_file: str, seed_id: str):
    command = "/Applications/MEGAN_/tools/gc-assembler --input " + input_file + \
              " -o " + output_file + \
              " -fun SEED -id " + seed_id + " -mic 60"
    os.system(command)
    print(command)


# PROTEINS THAT THE INPUT FILE IS SEARCHED FOR
# RIGHT NOW THESE ARE IMPORTANT ENZYMES INVOLVED IN THE REVERSE BETA OXIDASE (RBO) PATHWAY
proteins_of_interest = []

proteins_of_interest.append("butyry"
                            "l-CoA")
proteins_of_interest.append("acyl-CoA dehydrogenase")
proteins_of_interest.append("thiolase")
proteins_of_interest.append("thl")
proteins_of_interest.append("hydroxybutyrate")
proteins_of_interest.append("butyryl-CoA:acetate")

# CREATE DICT WITH CONTIGS AS KEYS AND BEST MATCHING PROTEIN AS VALUE

with open(INPUT_FILE) as file:
    for line in file.readlines():
        if ("Contig-" in line):
            # print(line.split()[1])
            contig_id = line.split()[1]
            # print(line)
        if ("WP_" in line and not line.startswith(">")):
            # print(line)
            contig_hit = line
        CONTIG_DICT[contig_id] = contig_hit

# ITERATE OVER ALL DICTIONARY VALUES (BLAST RESULT LINES CONTAINING PROTEINS AND IDENTITIES) AND PRINT PROTEINS OF CHOICE
for protein_count, protein in enumerate(CONTIG_DICT.values()):
    if protein is "":
        continue
    full_line = protein.split()
    # print(full_line)
    # EXTRACT PROTEIN NAME FROM RESULT LINE
    for count, word in enumerate(full_line):
        first_name_word_index = count
        if word.isdigit():
            last_name_word_index = count - 1
    name_list = full_line[1:last_name_word_index]
    full_name = ""
    for word in name_list:
        full_name += word + " "
        # EXTRACT IDENTITY PERCENTAGE FROM RESULT LINE
    identity = full_line[-1]
    identity_int = int(identity[0:-1])
    # ONLY PRINT PROTEINS WITH WANTED IDENTITY
    if identity_int > 95:
        for interesting_protein in proteins_of_interest:
            if interesting_protein in full_name:
                print("Found protein " + full_name + " for contig " + str(protein_count) + " with identity " + identity)
