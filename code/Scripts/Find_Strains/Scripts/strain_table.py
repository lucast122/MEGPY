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
import re
from difflib import SequenceMatcher

organisms = []
organism_dict = {}
enzymes = []


# simple method to find similar strings used to group similar names
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def make_strain_csv(strain_file):
    print("Computing strain spreadsheet for file: " + strain_file)
    # open file and find all organism names with regex
    with open(strain_file, "r") as file:
        for line in file.readlines():
            organism = re.findall('\[(.*?)\]', line)
            enzyme = re.findall(r'\|\s(.*?) \[', line)
            # remove duplicates
            enzyme = list(set(enzyme))
            organism = list(set(organism))

            # print(enzyme)
            for e in enzyme:
                # Filter hypothetical proteins
                if "hypothetical" not in e:
                    organism_dict[e.lower()] = organism

    # for count,e in enumerate(enzymes):
    #     similarity = similar(e,enzymes[count-1])
    #     if similarity > 0.8:

    with open(strain_file + ".csv", 'w') as out_file:
        for key in organism_dict.keys():
            out_file.write(key + ';')
        out_file.write('\n')
        for idx, item in enumerate(organism_dict.items()):
            for key in organism_dict.keys():
                try:
                    strain = organism_dict[key][idx]
                    out_file.write(strain)
                    out_file.write(';')
                except IndexError:
                    out_file.write(';')

            out_file.write('\n')


# compute strain csv for all strain files computed by previous scripts
dirs = [x[0] for x in os.walk(".")]

for d in dirs:
    files = os.listdir(d)
    for f in files:
        if f.startswith("strains") and not f.endswith(".csv") and not f.endswith(".py"):
            f = (d + "/" + f)
            make_strain_csv(f)
