import argparse

parser = argparse.ArgumentParser(description = "Calculates specificity scores for genomic regions with multiple tissue scores.")
parser.add_argument("input", type = str, help = "Path to the text file that lists the genomic regions on which tissue specificity will be calculated. \
                    Assumed format for the file: CHROMOSOME  START_COORDINATE    END_COORDINATE  [ADDITIONAL_COLUMNS]  TISSUE_1_SCORE  [TISSUE_2_SCORE] ... \
                    A header describing the columns is assumed.")
parser.add_argument("tissue", type = str, help = "Tissue to calculate specificity for.")
parser.add_argument("--ns", type = str, default = "_", help = "The separator symbol for extracting the tissue name from the column name (naming symbol). '_' by default. \
                    Example naming format for score columns: 'HL_SCORE' (the part after the separator must always be 'SCORE').")
parser.add_argument("--output", type = str, default = "stdout", help = "Path to the output file (text). Output will be in .BED format and preserve existing columns while adding new ones. If specified as not a proper name but 'stdout', which is the default option, the output will simply be printed and can be piped. \
                    Specificity score is defined as the average of differences between the score for the tissue in question and that for all others. \
                    It can take values between -1 and 1, and higher scores indicate higher specificity.")
args = parser.parse_args()

try:
    with open(args.input, "r") as input:
        input_lines = input.readlines()
except FileNotFoundError:
    print("The input file doesn't exist.")
    exit()

region_dict = []
header = []
dict_template = []

for i, line in enumerate(input_lines):
    if line[0] != "#":
        elts = line.strip().split()
        
        coords = [int(elts[1]), int(elts[2]), elts[0]]
        dict_construct = {}

        region_dict.append([coords, {}])
        for col in dict_template:
            dict_construct[col] = ""

            region_dict[i - 1][1] = dict_construct
            
            keys_vals = zip(list(dict_construct.keys()), elts[3:])

            for key, val in keys_vals:
                try:
                    val = float(val)
                except ValueError:
                    pass
                region_dict[i - 1][1][key] = val
        else:
            pass
    else:
        header = line.strip().split()
        
        for i, col in enumerate(header):
            if i != 0 and i != 1 and i != 2 and i != 3:
                dict_template.append(col)
            else:
                pass
        
        header.append("SPECIFICITY_SCORE")

main_key = args.tissue + args.ns + "SCORE"

for i, entry in enumerate(region_dict):
    differences = []
    average_diff = 0

    for key in entry[1].keys():
        if ("SCORE" in key) and (key != main_key):
            differences.append(float(entry[1][main_key]) - float(entry[1][key]))
        else:
            pass
    
    average_diff = sum(differences) / len(differences)

    region_dict[i][1]["SPECIFICITY_SCORE"] = average_diff

if args.output == "stdout":
    for i, col in enumerate(header):
        if col == "#":
            print(col, end = "")
        elif i == 1:
            print(" " + col, end = "")
        elif i == (len(header) - 1):
            print("\t" + col)
        else:
            print("\t" + col, end = "")
    
    for entry in region_dict:
        print(entry[0][2] + "\t" + str(int(entry[0][0])) + "\t" + str(int(entry[0][1])), end = "")

        for i, key in enumerate(list(entry[1].keys())):
            if (entry[1][key] == 0) or (entry[1][key] == 1):
                entry[1][key] = int(entry[1][key])
            else:
                pass

            if i == (len(list(entry[1].keys())) - 1):
                print("\t" + str(entry[1][key]))
            else:
                print("\t" + str(entry[1][key]), end = "")
else:
    with open(args.output, "w") as output:
        for i, col in enumerate(header):
            if col == "#":
                output.write(col)
            elif i == 1:
                output.write(" " + col)
            elif i == (len(header) - 1):
                output.write("\t" + col + "\n")
            else:
                output.write("\t" + col)

        for entry in region_dict:
            output.write(entry[0][2] + "\t" + str(int(entry[0][0])) + "\t" + str(int(entry[0][1])))

            for i, key in enumerate(list(entry[1].keys())):
                if (entry[1][key] == 0) or (entry[1][key] == 1):
                    entry[1][key] = int(entry[1][key])
                else:
                    pass
                 
                if i == (len(list(entry[1].keys())) - 1):
                    output.write("\t" + str(entry[1][key]) + "\n")
                else:
                    output.write("\t" + str(entry[1][key]))
