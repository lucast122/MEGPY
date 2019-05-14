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

import os

PATHS = ["/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW104_TCCTGAGCGTAAGGAG_L008_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW105_GGACTCCTGTAAGGAG_L001_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW105_GGACTCCTGTAAGGAG_L008_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW106_TAGGCATGGTAAGGAG_L001_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW106_TAGGCATGGTAAGGAG_L008_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW107_CTCTCTACGTAAGGAG_L001_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW107_CTCTCTACGTAAGGAG_L008_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW108_CAGAGAGGGTAAGGAG_L001_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW108_CAGAGAGGGTAAGGAG_L008_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW109_GCTACGCTGTAAGGAG_L001_R1_001_trimmed..daa",
         "/Volumes/Elements/Uni/Data_Master_Thesis/Mapping/JW109_GCTACGCTGTAAGGAG_L008_R1_001_trimmed..daa"]

# run gc_assembly_map for all given daa files

for path in PATHS:
    command = "python3 gc_assembly_map.py --input " + path + " --ids 441,42,56"
    print(command)
    os.system(command)
