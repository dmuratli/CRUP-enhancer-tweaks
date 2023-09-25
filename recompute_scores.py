import argparse
import sys

def find_coord(coords, chr, entry_list, start_pos = 0):
    for i in range(start_pos, len(entry_list)):
        if (coords[0] <= entry_list[i][0] < coords[1]) and (coords[0] < entry_list[i][1] <= coords[1]) and (entry_list[i][3] == chr):
            return i
        else:
            pass
        
    return sys.maxsize

parser = argparse.ArgumentParser(description = "Recomputes region scores based on the inputted bins.")
parser.add_argument("bins", type = str, nargs="*", help = "Path(s) to the text file(s) that contain(s) all the bins along with their scores. \
                    .BED format is assumed with an extra column of scores/probabilities produced for the output. Can take one or more inputs - if multiple files are provided, \
                    multiple scores are calculated and fittingly described.")
parser.add_argument("regions", type = str, help = "Path to the text file that lists predicted regions. \
                    The file must have its previous score columns spliced out beforehand to prevent column duplications. \
                    .BED format is assumed (extra columns are allowed). Using a header to describe the columns is \
                    a requirement when the file has columns not included by default in the .BED format. \
                    Can only take one input.")
parser.add_argument("--output", type = str, default = "stdout", help = "Path to the output file (text). Output will be in .BED format and preserve or potentially add extra columns. If specified as not a proper name but 'stdout', which is the default option, the output will simply be printed and can be piped.")
parser.add_argument("--ns", type = str, default = "_", help = "The separator symbol for extracting the tissue name from the file name (naming symbol). '_' by default.")
args = parser.parse_args()

region_dict = {}
header = []

try:
    with open(args.regions, "r") as regions:
        regions_lines = regions.readlines()

        for line in regions_lines:
            if line[0] != "#":
                elts = line.strip().split()
                coords = (int(elts[1]), int(elts[2]), elts[0])

                if coords not in region_dict:
                    region_dict[coords] = []
                    
                    for elt in elts[3:]: # None of the columns here should contain the previous score for the region!
                        region_dict[coords].append(elt)
                else:
                    pass
            else:
                header = line.strip().split()
except FileNotFoundError:
    print("The region file does not exist.")
    exit()

if len(header) == 0:
    header = ["#", "CHROMOSOME", "START_COORDINATE", "END_COORDINATE"]
else:
    pass

for i in range(0, len(args.bins)):
    try:
        with open(args.bins[i], "r") as bins:
            bins_lines = bins.readlines()
            bins_list = []

            for line in bins_lines:
                if line[0] != "#":
                    elts = line.strip().split()
                    ord_lst = [int(elts[1]), int(elts[2]), float(elts[3]), elts[0]]
                    bins_list.append(ord_lst)
                else:
                    pass
    except FileNotFoundError:
        print("The bin file does not exist.")
        exit()

    for region in list(region_dict.keys()):
        coords = [region[0], region[1]]
        chr = region[2]
        probs = []

        pos = find_coord(coords, chr, bins_list)

        if pos != sys.maxsize:
            while (coords[0] <= bins_list[pos][0] < bins_list[pos][1] <= coords[1]):
                probs.append(bins_list[pos][2])
                pos += 1
        else:
            probs.append(0)

        highest_prob = max(probs)

        region_dict[region].append(highest_prob)

if args.output == "stdout":
    tissue_names = []

    for bin_file in args.bins:
        elts_1 = bin_file.split("/")
        elts_2 = elts_1[-1].split(args.ns)
        tissue_names.append(elts_2[0])

    for tissue_name in tissue_names:
        modified_tissue_name = tissue_name + "_SCORE"
        header.append(modified_tissue_name)

    for i, col in enumerate(header):
        if i == 0:
            print(col, end = "")
        elif i == 1:
            print(" " + col, end = "")
        elif i == (len(header) - 1):
            print("\t" + col)
        else:
            print("\t" + col, end = "")

    for region in list(region_dict.keys()):
        print(region[2] + "\t" + str(region[0]) + "\t" + str(region[1]), end = "")

        for i, val in enumerate(region_dict[region]):
            if i != (len(region_dict[region]) - 1):
                print("\t" + str(val), end = "")
            else:
                print("\t" + str(val))
else:
    with open(args.output, "w") as updated:
        tissue_names = []

        for bin_file in args.bins:
            elts_1 = bin_file.split("/")
            elts_2 = elts_1[-1].split(args.ns)
            tissue_names.append(elts_2[0])

        for tissue_name in tissue_names:
            modified_tissue_name = tissue_name + "_SCORE"
            header.append(modified_tissue_name)

        for i, col in enumerate(header):
            if i == 0:
                updated.write(col)
            elif i == 1:
                updated.write(" " + col)
            elif i == (len(header) - 1):
                updated.write("\t" + col + "\n")
            else:
                updated.write("\t" + col)

        for region in list(region_dict.keys()):
            updated.write(region[2] + "\t" + str(region[0]) + "\t" + str(region[1]))

            for i, val in enumerate(region_dict[region]):
                if i != (len(region_dict[region]) - 1):
                    updated.write("\t" + str(val))
                else:
                    updated.write("\t" + str(val) + "\n")
