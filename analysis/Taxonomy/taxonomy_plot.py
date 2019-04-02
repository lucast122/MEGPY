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
"""

import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('threshold', type=float,
                    help='Float vlue for plotting cutoff. Only classes with higher percentage than this value are shown')

parser.add_argument('-p', '--plot', action='store_true',
                    help='If set MEGAN is not run and just the plot is created')
args = parser.parse_args()

# Threshold in percent for plotting. Every taxonomic class with a percentage higher than this is used for the plot
threshold = args.threshold
mapping_file = ""
mapping_file_path = "/Volumes/Elements/Uni/Data Master Thesis/ready/"
count_folder_path = "/Users/Timo/Dropbox/Thesis/METAMAP/analysis/Count Data/"
count_file_name = ""

phyla = []
phyla_combined = []
values = []
pos_count = 0
pos = []
data_bars = []

# Create commands.txt with the proper commands for the right file

os.chdir(mapping_file_path)
# For loop over all .daa files in folder specified by mapping_file_path
# Creates a specific command.txt for each .daa file present
for file in glob.glob("*.daa"):
    mapping_file = file[:-4]
    print("Current mapping file: " + str(mapping_file) + ".daa")
    os.chdir(count_folder_path)
    with open(str(mapping_file) + '_commands.txt', 'w') as the_file:
        count_file_name = str(mapping_file) + "_taxon_to_percent.txt"
        count_file_absolute_path = count_folder_path + count_file_name
        the_file.write("open viewer=Taxonomy;\n")
        the_file.write("open file='" +
                       str(mapping_file_path) + str(mapping_file) + ".daa';\n")

        the_file.write("uncollapse nodes = all;\n")

        the_file.write("select rank=Order;\n")
        the_file.write("select nodes=subtree;\n")
        the_file.write("export what=CSV format=taxonName_to_percent separator=tab counts=assigned file='" +
                       str(count_file_absolute_path) + "';\n")
        the_file.write("quit;")

command = "/Applications/MEGAN_/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E < " + \
          str(mapping_file) + "_commands.txt"
if (not args.plot):
    os.system(command)

# find all taxon_to_percent  .txt
# files in the taxon_to_percent folder
# and use data for plotting. For loop to iterate over all files
os.chdir(count_folder_path)

for file in glob.glob("*taxon_*.txt"):
    phyla = []
    values = []
    phyla_with_percentage = []

    file_name = file

    file_ID = os.path.splitext(file_name)[0]
    full_path = count_folder_path + file_name

    print("Creating taxonomic plot for file " + file_ID)

    with open(full_path, "r") as f:
        content = f
        array = []
        for line in f:
            array.append(line)

    for entry in array:
        entry = str(entry)
        phylum = entry.split("\t")[0]
        value = float(entry.split("\t")[1].rstrip())
        if value > threshold:
            phyla.append(phylum[1:-1])
            values.append(value)

    a = 0
    for value in values:
        a += value

    other = 100 - a
    values.append(other)
    # round values to 2 decimal numbers
    values = ['%.2f' % elem for elem in values]

    c = 0
    for phy in phyla:
        phy = phy + " [" + str(values[c]) + "%]"
        phyla_with_percentage.append(phy)
        c += 1

    # add category to plot that represents everything below threshold

    phyla_with_percentage.append("Other " + "[" + str(values[-1]) + "%]")
    phyla.append("Other")
    phyla_combined += phyla
    data = values
    data_bars += values
    recipe = phyla_with_percentage
    fig, ax = plt.subplots(figsize=(16, 10), subplot_kw=dict(aspect="equal"))

    wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=-40)

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)

    kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(recipe[i], xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw)

    # Title messses picture up. Use latex fig title instead
    # ax.set_title("Distribution of taxa on phylum level")
    try:
        fig.savefig(count_folder_path + file_ID + ".pdf", bbox_inches='tight')
        print("Plot successfully created for " + str(file_ID))
        plt.close()
        pass
    except Exception as e:
        print("Could not create plot successfully.")
        plt.close()
    else:
        pass
    finally:
        pass

# start 2nd plot

print("Creating grouped bar chart ...")
bar_width = 0.005
# plt.bar(pos, data, bar_width, color='blue', edgecolor='black')
for dat in data_bars:
    pos.append(pos_count)
    pos_count += 0.0055
data_int = [round(int(float(i))) for i in data_bars]

count = 0

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

for dat in data_bars:
    bar = plt.bar(pos, height=data_int, width=bar_width, label=phyla_combined[count],color=colors)
    count += 1

plt.xticks()
plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
plt.margins(x=0.03)
plt.xlim(0)
plt.xlabel('Taxonomic Class', fontsize=12)
plt.ylabel('Relative abundance in percent', fontsize=12)
plt.title('Taxonomic analysis of the dataset ' + str(mapping_file), fontsize=12)

plot_name = str(mapping_file)+"_grouped_barchart.pdf"

try:
    plt.savefig(plot_name, bbox_inches="tight")
    print("Saved plot as " + plot_name)
    plt.close()
    pass
except Exception as e:
    print("Could not create plot successfully.")
    plt.close()
else:
    pass
finally:
    pass


plt.show()

