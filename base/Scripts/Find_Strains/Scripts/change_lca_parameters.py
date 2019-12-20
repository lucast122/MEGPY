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

import argparse
import glob
import ntpath
import os

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('-i', '--input', type=str, required=True,
                    help='Path to input directory containing the mapping (.daa) files')
parser.add_argument('-o', '--out', type=str, required=True,
                    help='Path to output directory where intermediate results are written to')

args = parser.parse_args()

# change MEGAN_COMMAND string to fit your MEGAN path
MEGAN_COMMAND = "/Applications/MEGAN_/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E < "


def change_min_identity():
    mapping_file_path = args.input

    for file in glob.glob(mapping_file_path + "*.daa"):
        mapping_file_name = ntpath.basename(file)
        # mapping_file_name = "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW101_L001_R0_001.daa"
        print("Current mapping file: " + mapping_file_name)
        if not os.path.isdir(args.out):
            os.mkdir(args.out)
        with open(args.out + mapping_file_name + '_commands.txt', 'w') as the_file:
            count_file_name = str(mapping_file_name) + "_taxon_to_percent.txt"
            the_file.write("open viewer=Taxonomy;\n")
            the_file.write("open file='" +
                           str(file) + "';\n")

            the_file.write("recompute minSupportPercent = 0.01 minSupport = 1 minScore = 50.0 maxExpected = 0.01"
                           " minPercentIdentity = 80.0 topPercent = 10.0 lcaAlgorithm = naive"
                           " lcaCoveragePercent = 100.0 minPercentReadToCover = 0.0"
                           " minPercentReferenceToCover = 0.0 minComplexity = 0.0 longReads = false"
                           " pairedReads = false useIdentityFilter = false readAssignmentMode"
                           " = readCount fNames = SEED;\n")

            the_file.write("quit;")
            the_file.close()

            command = MEGAN_COMMAND + \
                      args.out + str(mapping_file_name) + "_commands.txt"

            os.system(command)


def main():
    change_min_identity()


if __name__ == '__main__':
    main()
