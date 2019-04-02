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


This script is the second part of the analysis pipeline for finding specific strains in a metagenomics dataset that are
associated with a specific gene.

After successfully running gc_assembly_map.py and getting a gene-centric assembly and blast output this script takes a
list of strains of interest and searches the blast result file for those strains. It then returns a detailed description
about the strains present in the dataset.
"""

import argparse

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('strains', type=str,
                    help='File containing strains of interest. Use a plain text file and use a new line for each strain')
parser.add_argument('blast_results', type=str,
                    help='File with the blast results created with gc_assembly_map.py')
args = parser.parse_args()
