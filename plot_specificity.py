import seaborn as sns
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = "Calculates specificity scores for genomic regions with multiple tissue scores.")
parser.add_argument("input", type = str, help = "Path to the text file that lists the genomic regions with their respective specificity scores. \
                    Assumed format for the file: CHROMOSOME  START_COORDINATE    END_COORDINATE ACTIVE_TISSUES [ADDITIONAL_COLUMNS]  SPECIFICITY_SCORE ... \
                    A header describing the columns is assumed.")
parser.add_argument("--tissues", type = str, default = "ALL", help = "Name(s) of the tissue(s) whose associated specificity scores will be plotted. \
                    If multiple tissues are going to be provided, they have to be separated by a comma (e.g. 'FL,HL,B'). \
                    By default, no filtering is done - specificity scores of regions labelled with any or all of the tissues are included in the plot.")
parser.add_argument("--sep", type = str, default = ",", help = "The separator for tissue names under 'ACTIVE_TISSUES'. ',' by default.")
args = parser.parse_args()

try:
    with open(args.input, "r") as input:
        input_lines = input.readlines()
except FileNotFoundError:
    print("The input file doesn't exist.")
    exit()

scores = []
temp_dict = {}
header = []

for line in input_lines:
    if line[0] != "#":
        elts = line.strip().split()

        for i, elt in enumerate(elts[3:]):
            try:
                val = float(elt)
            except ValueError:
                val = elt
        
            temp_dict[list(temp_dict.keys())[i]] = val
        
        if args.tissues != "ALL":
            if sorted(temp_dict["ACTIVE_TISSUES"].strip().split(",")) == sorted((args.tissues).strip().split(",")):
                scores.append(temp_dict["SPECIFICITY_SCORE"])
            else:
                pass
        else:
            scores.append(temp_dict["SPECIFICITY_SCORE"])
    else:
        header = line.strip().split()

        for i, col in enumerate(header):
            if i != 0 and i != 1 and i != 2 and i != 3:
                temp_dict[col] = 0
            else:
                pass

sns.set_style("darkgrid")

plot1 = sns.histplot(scores)
plt.savefig(args.input + "_histplot.jpg", dpi = 300)
plt.close()

plot2 = sns.kdeplot(scores)
plt.savefig(args.input + "_kdeplot.jpg", dpi = 300)
plt.close()