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
import random
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Parse command line arguments.')
parser.add_argument('-t', '--threshold', type=float,
                    help='Float value for plotting cutoff. Only classes with higher '
                         'percentage than this value are shown')

parser.add_argument('-e', '--extract', action='store_true',
                    help='Provide this parameter to extract the taxon_to_percent data from the provided .daa files using MEGAN'
                         'for this to work MEGAN ultimate has to be installed and ')

parser.add_argument('--donut', action='store_true',
                    help='Provide this parameter to create the donut plots ')

parser.add_argument('--bar', action='store_true',
                    help='Provide this parameter to create the grouped bar plots')

parser.add_argument('-d', '--daa', type=str,
                    help='Path to the folder containing the .daa files of interest')

parser.add_argument('-o', '--output', type=str,
                    help='Path to output folder')

args = parser.parse_args()

"""
CUSTOMIZE MEGAN_COMMAND:
YOU MUST CHANGE THE PATH TO MEGAN ULTIMATE TO FIT YOUR SYSTEM
"""
MEGAN_COMMAND = "/Applications/MEGAN_/MEGAN.app/Contents/MacOS/JavaApplicationStub -g -E < "

# Threshold in percent for plotting. Every taxonomic class with a percentage higher than this is used for the plot
threshold = args.threshold
mapping_file = ""
mapping_file_path = args.daa
count_folder_path = ""
count_file_name = ""

out_path = args.output
daa_names_trimmed = []
phyla = []
phyla_combined = []
values = []
pos = []
temp_bars = []
data_bars = []
phyla_pairwise = []
data_pairwise = []
temp = []
xticks_pos = [0]
temp2 = []
data_donut_plots = []
labels_donut_plots = []
daa_file_names = []
# Create commands.txt with the proper commands for the right file

# For loop over all .daa files in folder specified by mapping_file_path
# Creates a specific command.txt for each .daa file present,

# find all taxon_to_percent  .txt
# files in the taxon_to_percent folder
# and use data for plotting. For loop to iterate over all files


for count, file in enumerate(glob.glob("Output/" + "*taxon_*.txt"), 1):
    daa_file_names.append(file)
    phyla = []
    values = []
    phyla_with_percentage = []

    file_name = file

    file_ID = os.path.splitext(file_name)[0]
    full_path = count_folder_path + file_name

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
    temp_bars += values

    data_bars.append(0)

    phyla_pairwise.append(phyla)
    data_pairwise.append(values)

    data_donut_plots.append(data)
    labels_donut_plots.append(phyla_with_percentage)


# call important methods


def main():
    if args.extract:
        extract_taxon_to_percent()
    if args.donut:
        for count, lab in enumerate(labels_donut_plots):
            name = create_donut_plots(data_donut_plots[count], lab, daa_file_names[count])
            daa_names_trimmed.append(name)
    if args.bar:
        create_grouped_barplot(data_bars, daa_names_trimmed)


def extract_taxon_to_percent():
    for file in glob.glob(mapping_file_path + "*.daa"):
        mapping_file_name = ntpath.basename(file)
        print("Current mapping file: " + mapping_file_name)

        with open("Output/" + mapping_file_name + '_commands.txt', 'w') as the_file:
            count_file_name = str(mapping_file_name) + "_taxon_to_percent.txt"
            count_file_absolute_path = count_folder_path + count_file_name
            the_file.write("open viewer=Taxonomy;\n")
            the_file.write("open file='" +
                           str(file) + "';\n")

            the_file.write("uncollapse nodes = all;\n")

            the_file.write("select rank=Order;\n")
            the_file.write("select nodes=subtree;\n")
            the_file.write(
                "export what=CSV format=taxonName_to_percent separator=tab counts=assigned file='" + out_path +
                str(count_file_absolute_path) + "';\n")
            the_file.write("quit;")
            the_file.close()

            command = MEGAN_COMMAND + \
                      "Output/" + str(mapping_file_name) + "_commands.txt"
        if not args.plot:
            try:
                os.system(command)
            except OSError:
                print("Failed to start MEGAN. Please configure the correct path")


def create_donut_plots(data_donuts, labels, name):
    name = name.split('/')[1].split("_")[0] + name.split('/')[1].split("_")[1] + name.split('/')[1].split("_")[2] + \
           name.split('/')[1].split("_")[3]

    my_colors = [colors(n)[0] for n, _ in enumerate(labels, 1)]
    color_dict = dict(zip(labels, my_colors))

    recipe = labels
    fig, ax = plt.subplots(figsize=(16, 10), subplot_kw=dict(aspect="equal"))
    print("Creating taxonomic plot for " + name)
    wedges, texts = ax.pie(data_donuts, wedgeprops=dict(width=0.5), startangle=-40, colors=my_colors)
    #
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    #
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
        fig.savefig(name + "_taxon_donut" ".pdf", bbox_inches='tight')
        print("Plot successfully created for " + name)
        plt.close()
        pass
    except Exception as e:
        print("Could not create plot successfully.")
        plt.close()
    else:
        pass
    finally:
        pass

    return (name)


def create_grouped_barplot(data, daa_names):
    # sns.set_palette(sns.color_palette("hls", 20))
    pos_count = 3
    print("Creating grouped bar chart ...")
    bar_width = 3

    for dat in data:
        if dat == 0:
            pos_count += 20
            xticks_pos.append(pos_count)
            continue
        pos.append(pos_count)
        pos_count += 6

    dat = temp_bars

    data_int = [round(int(float(i))) for i in dat]

    keys = range(len(data_int))
    phyla_dict = dict(zip(keys, phyla_combined))
    phyla_non_redundant = list(dict.fromkeys(phyla_combined))
    # create dictionary mapping taxon labels to colors
    my_colors = [colors(n)[0] for n, _ in enumerate(phyla_non_redundant, 1)]
    color_dict = dict(zip(phyla_non_redundant, my_colors))
    bar_container = plt.bar(pos, height=data_int, width=bar_width)
    # set bar color according to taxon class

    labels = []
    for counter, b in enumerate(bar_container):
        lab = phyla_dict[counter]
        labels.append(lab)
        b.set_label(lab)
        col = color_dict[lab]
        b.set_color(col)
    # labels = list(dict.fromkeys(labels))
    # remove duplicates from legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))

    # create legend
    plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04, 1), loc="upper left")
    # plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.xticks(xticks_pos, daa_names, rotation="vertical")
    plt.margins(x=0.03)

    plt.xlim(1)
    plt.xlabel('Dataset ID', fontsize=12)
    plt.ylabel('Relative abundance in percent', fontsize=12)
    plt.title('Taxonomic analysis of the dataset ' + str(mapping_file), fontsize=12)

    plot_name = "grouped_barchart.pdf"
    # plt.show()
    print("Saving plot ...")
    plt.savefig(plot_name, bbox_inches="tight")
    print("Saved plot as " + plot_name)
    plt.close()


def colors(n):
    if n == 0:
        return
    ret = []
    r = int(random.random() * 256)
    g = int(random.random() * 256)
    b = int(random.random() * 256)
    assert n != 0
    step = 256 / n
    for i in range(n):
        r += step
        g += step
        b += step
        r = (int(r) % 256) / 256
        g = (int(g) % 256) / 256
        b = (int(b) % 256) / 256
        ret.append((r, g, b))
        return ret


if __name__ == '__main__':
    main()
