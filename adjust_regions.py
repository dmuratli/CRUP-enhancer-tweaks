import sys
import argparse

parser = argparse.ArgumentParser(description = "Adjusts regions by shrinking and/or extending them based on their underlying scores and user-configured cut-offs.")
parser.add_argument("bins", type = str, help = "Path to the text file that contains all the bins along with their scores. \
                    Assumed format for the file: CHROMOSOME START_COORDINATE   END_COORDINATE  SCORE/PROBABILITY")
parser.add_argument("regions", type = str, help = "Path to the text file that lists predicted regions. \
                    Assumed format for the file is the same as that of bins.")
parser.add_argument("output", type = str, help = "Path to the output file (text). \
                    Output format will be the same as that of bins and regions in the input.")
parser.add_argument("--operation", type = str, default = "both", help = "Operation(s) that will be performed on the regions. Options are 'both', 'shrinkage', 'extension'. Set to 'both' by default.")
parser.add_argument("--extension_cutoff", type = str, default = "0.5", help = "Cut-off value for extending regions - regions will be extended if adjacent bins have scores higher than or equal to this value. \
                    Defined as 0.5 by default. You can either set an absolute cut-off (e.g. 0.5) or dynamically define it in terms of the highest score within the region by specifying operations and values \
                    (e.g. /n*x, where n and x are numerical values; this would divide the highest score by n and multiply the result by x). Parentheses and spaces between characters are NOT allowed in the expression for setting a dynamic cut-off.")
parser.add_argument("--shrinkage_cutoff", type = str, default = "/2", help = "Cut-off value for shrinking regions - regions will be shortened as long as there are bins that have strictly lower scores than this value on the edges. \
                    Dynamically defined as 'highest score within the region / 2' by default. Similar to the extension cut-off, it can be defined as an absolute number or relative to the highest score exactly the same way described above.")
parser.add_argument("--margin", type = float, default = 0, help = "Minimum length allowed for regions. 'Separate' regions will be merged with each other if the distance between them is smaller than or equal to the margin. \
                    Furthermore, if a region is smaller than or equal to the margin in size but have no adjacent regions, it will simply be removed. \
                    0 (turned off) by default.")
args = parser.parse_args()

def find_coord(direction, operation, coord, chr, entry_list, start_pos = 0):
    if operation == "extend":
        for i in range(start_pos, len(entry_list)):
            if (direction == "left") and (entry_list[i][0] < coord <= entry_list[i][1]) and (chr == entry_list[i][3]):
                return i
            elif (direction == "right") and (entry_list[i][0] <= coord < entry_list[i][1]) and (chr == entry_list[i][3]):
                return i
            else:
                pass
    elif operation == "shrink":
        for i in range(start_pos, len(entry_list)):
            if (direction == "left") and (entry_list[i][0] <= coord < entry_list[i][1]) and (chr == entry_list[i][3]):
                return i
            elif (direction == "right") and (entry_list[i][0] < coord <= entry_list[i][1]) and (chr == entry_list[i][3]):
                return i
            else:
                pass
    else:
        pass
        
    return sys.maxsize

def divide(s, val): # Adapted from https://stackoverflow.com/questions/43389684/how-can-i-split-a-string-of-a-mathematical-expressions-in-python
    r = [""]
    s = str(val) + s

    for i in s.replace(" ", ""):
        if i.isdigit() and r[-1].isdigit():
            r[-1] = r[-1] + i
        else:
            r.append(i)
    
    return r[1:]

def interpret_cutoff(input_cutoff, highest_prob):
    try:
        cutoff = float(input_cutoff)
        
        return cutoff
    except:
        cutoff = input_cutoff

        exp = divide(cutoff, highest_prob)

        running_val = 0
        prev_op = ""

        for i in range(0, len(exp)):
            if exp[i].isdigit():
                if i == 0:
                    running_val = float(exp[i])
                else:
                    if prev_op == "+":
                        running_val += float(exp[i])
                    elif prev_op == "-":
                        running_val -= float(exp[i])
                    elif prev_op == "*":
                        running_val *= float(exp[i])
                    elif prev_op == "/":
                        running_val /= float(exp[i])
            else:
                prev_op = exp[i]
        
        return running_val

try:
    with open(args.bins, "r") as bins:
        bins_list = []
        bins_lines = bins.readlines()

        for line in bins_lines:
            if line[0] != "#":
                elts = line.strip().split()
                entry = [int(elts[1]), int(elts[2]), float(elts[3]), elts[0]]
                bins_list.append(entry)
            else:
                pass
except FileNotFoundError:
    print("The bin file does not exist.")
    exit()

try:
    regions = open(args.regions, "r")
except FileNotFoundError:
    print("The region file does not exist.")
    exit()

refined = open(args.output, "w")

refined.write("# CHROMOSOME START_COORDINATE   END_COORDINATE  SCORE/PROBABILITY\n")

if args.operation == "both":
    for line in regions:
        if line[0] != "#":
            elts = line.strip().split()
            refined.write(str(elts[0]) + "\t")
            highest_prob = float(elts[3])
            new_elts = [int(elts[1]), int(elts[2]), highest_prob, elts[0]]

            shr_cutoff = interpret_cutoff(args.shrinkage_cutoff, highest_prob)
            ext_cutoff = interpret_cutoff(args.extension_cutoff, highest_prob)

            print("Modifying the region [" + str(new_elts[0]) + ", " + str(new_elts[1]) + ") on " + new_elts[3] + ":")

            shr_idx1 = find_coord("left", "shrink", new_elts[0], elts[0], bins_list)

            if shr_idx1 != sys.maxsize:
                shr_idx2 = find_coord("right", "shrink", new_elts[1], elts[0], bins_list, shr_idx1)
            else:
                shr_idx2 = find_coord("right", "shrink", new_elts[1], elts[0], bins_list)
            
            while (shr_idx1 != sys.maxsize) or (shr_idx2 != sys.maxsize):
                if shr_idx1 != sys.maxsize:
                    print("Shrinking from the left...")
                    while (bins_list[shr_idx1][0] <= new_elts[0] < bins_list[shr_idx1][1]) and (bins_list[shr_idx1][2] < shr_cutoff) and (bins_list[shr_idx1][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[shr_idx1][0]) + ", " + str(bins_list[shr_idx1][1]) + ")...")
                        new_elts[0] = bins_list[shr_idx1][1]

                        shr_idx1 += 1
                    
                    shr_idx1 = sys.maxsize
                else:
                    pass

                if shr_idx2 != sys.maxsize:
                    print("Shrinking from the right...")
                    while (bins_list[shr_idx2][0] < new_elts[1] <= bins_list[shr_idx2][1]) and (bins_list[shr_idx2][2] < shr_cutoff) and (bins_list[shr_idx2][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[shr_idx2][0]) + ", " + str(bins_list[shr_idx2][1]) + ")...")
                        new_elts[1] = bins_list[shr_idx2][0]

                        shr_idx2 -= 1
                    
                    shr_idx2 = sys.maxsize
                else:
                    pass

            ext_idx1 = find_coord("left", "extend", new_elts[0], elts[0], bins_list)

            if ext_idx1 != sys.maxsize:
                ext_idx2 = find_coord("right", "extend", new_elts[1], elts[0], bins_list, ext_idx1)
            else:
                ext_idx2 = find_coord("right", "extend", new_elts[1], elts[0], bins_list)

            while (ext_idx1 != sys.maxsize) or (ext_idx2 != sys.maxsize):
                if ext_idx1 != sys.maxsize:
                    print("Extending to the left...")
                    while (bins_list[ext_idx1][0] < new_elts[0] <= bins_list[ext_idx1][1]) and (bins_list[ext_idx1][2] >= ext_cutoff) and (bins_list[ext_idx1][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[ext_idx1][0]) + ", " + str(bins_list[ext_idx1][1]) + ")...")
                        new_elts[0] = bins_list[ext_idx1][0]

                        if bins_list[ext_idx1][2] > highest_prob:
                            highest_prob = bins_list[ext_idx1][2]
                        else:
                            pass

                        ext_idx1 -= 1
                    
                    ext_idx1 = sys.maxsize
                else:
                    pass

                if ext_idx2 != sys.maxsize:
                    print("Extending to the right...")
                    while (bins_list[ext_idx2][0] <= new_elts[1] < bins_list[ext_idx2][1]) and (bins_list[ext_idx2][2] >= ext_cutoff) and (bins_list[ext_idx2][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[ext_idx2][0]) + ", " + str(bins_list[ext_idx2][1]) + ")...")
                        new_elts[1] = bins_list[ext_idx2][1]
                        
                        if bins_list[ext_idx2][2] > highest_prob:
                            highest_prob = bins_list[ext_idx2][2]
                        else:
                            pass

                        ext_idx2 += 1
                    
                    ext_idx2 = sys.maxsize
                else:
                    pass
                
            if (highest_prob == 0) or (highest_prob == 1):
                highest_prob = int(highest_prob)
            else:
                pass

            refined.write(str(new_elts[0]) + "\t" + str(new_elts[1]) + "\t" + str(highest_prob) + "\n")
        else:
            pass
elif args.operation == "shrinkage":
    for line in regions:
        if line[0] != "#":
            elts = line.strip().split()
            refined.write(str(elts[0]) + "\t")
            highest_prob = float(elts[3])
            new_elts = [int(elts[1]), int(elts[2]), highest_prob, elts[0]]

            shr_cutoff = interpret_cutoff(args.shrinkage_cutoff, highest_prob)

            print("Modifying the region [" + str(new_elts[0]) + ", " + str(new_elts[1]) + ") on " + new_elts[3] + ":")

            shr_idx1 = find_coord("left", "shrink", new_elts[0], elts[0], bins_list)

            if shr_idx1 != sys.maxsize:
                shr_idx2 = find_coord("right", "shrink", new_elts[1], elts[0], bins_list, shr_idx1)
            else:
                shr_idx2 = find_coord("right", "shrink", new_elts[1], elts[0], bins_list)
            
            while (shr_idx1 != sys.maxsize) or (shr_idx2 != sys.maxsize):
                if shr_idx1 != sys.maxsize:
                    print("Shrinking from the left...")
                    while (bins_list[shr_idx1][0] <= new_elts[0] < bins_list[shr_idx1][1]) and (bins_list[shr_idx1][2] < shr_cutoff) and (bins_list[shr_idx1][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[shr_idx1][0]) + ", " + str(bins_list[shr_idx1][1]) + ")...")
                        new_elts[0] = bins_list[shr_idx1][1]

                        shr_idx1 += 1
                    
                    shr_idx1 = sys.maxsize
                else:
                    pass

                if shr_idx2 != sys.maxsize:
                    print("Shrinking from the right...")
                    while (bins_list[shr_idx2][0] < new_elts[1] <= bins_list[shr_idx2][1]) and (bins_list[shr_idx2][2] < shr_cutoff) and (bins_list[shr_idx2][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[shr_idx2][0]) + ", " + str(bins_list[shr_idx2][1]) + ")...")
                        new_elts[1] = bins_list[shr_idx2][0]

                        shr_idx2 -= 1
                    
                    shr_idx2 = sys.maxsize
                else:
                    pass

            if (highest_prob == 0) or (highest_prob == 1):
                highest_prob = int(highest_prob)
            else:
                pass

            refined.write(str(new_elts[0]) + "\t" + str(new_elts[1]) + "\t" + str(highest_prob) + "\n")
        else:
            pass
elif args.operation == "extension":
    for line in regions:
        if line[0] != "#":
            elts = line.strip().split()
            refined.write(str(elts[0]) + "\t")
            highest_prob = float(elts[3])
            new_elts = [int(elts[1]), int(elts[2]), highest_prob, elts[0]]

            ext_cutoff = interpret_cutoff(args.extension_cutoff, highest_prob)

            print("Modifying the region [" + str(new_elts[0]) + ", " + str(new_elts[1]) + ") on " + new_elts[3] + ":")

            ext_idx1 = find_coord("left", "extend", new_elts[0], elts[0], bins_list)

            if ext_idx1 != sys.maxsize:
                ext_idx2 = find_coord("right", "extend", new_elts[1], elts[0], bins_list, ext_idx1)
            else:
                ext_idx2 = find_coord("right", "extend", new_elts[1], elts[0], bins_list)

            while (ext_idx1 != sys.maxsize) or (ext_idx2 != sys.maxsize):
                if ext_idx1 != sys.maxsize:
                    print("Extending to the left...")
                    while (bins_list[ext_idx1][0] < new_elts[0] <= bins_list[ext_idx1][1]) and (bins_list[ext_idx1][2] >= ext_cutoff) and (bins_list[ext_idx1][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[ext_idx1][0]) + ", " + str(bins_list[ext_idx1][1]) + ")...")
                        new_elts[0] = bins_list[ext_idx1][0]

                        if bins_list[ext_idx1][2] > highest_prob:
                            highest_prob = bins_list[ext_idx1][2]
                        else:
                            pass

                        ext_idx1 -= 1
                    
                    ext_idx1 = sys.maxsize
                else:
                    pass

                if ext_idx2 != sys.maxsize:
                    print("Extending to the right...")
                    while (bins_list[ext_idx2][0] <= new_elts[1] < bins_list[ext_idx2][1]) and (bins_list[ext_idx2][2] >= ext_cutoff) and (bins_list[ext_idx2][3] == new_elts[3]):
                        print("Checking [" + str(bins_list[ext_idx2][0]) + ", " + str(bins_list[ext_idx2][1]) + ")...")
                        new_elts[1] = bins_list[ext_idx2][1]
                        
                        if bins_list[ext_idx2][2] > highest_prob:
                            highest_prob = bins_list[ext_idx2][2]
                        else:
                            pass

                        ext_idx2 += 1
                    
                    ext_idx2 = sys.maxsize
                else:
                    pass
                
            if (highest_prob == 0) or (highest_prob == 1):
                highest_prob = int(highest_prob)
            else:
                pass

            refined.write(str(new_elts[0]) + "\t" + str(new_elts[1]) + "\t" + str(highest_prob) + "\n")
        else:
            pass
else:
    print("Incorrect operation! Only options are 'shrinkage', 'extension', and 'both'.")
    exit()

refined.close()
regions.close()

if args.margin != 0:
    with open(args.output, "r") as refined:
        refined_lines = refined.readlines()
    
    lst_refined_lines = []

    for line in refined_lines:
        if line[0] != "#":
            line = line.strip().split()
            lst_refined_lines.append(line)
        else:
            pass

    with open(args.output, "w") as new_refined:
        new_refined.write("# CHROMOSOME START_COORDINATE   END_COORDINATE  SCORE/PROBABILITY\n")

        for i in range(0, len(lst_refined_lines)):
            c_chrom = lst_refined_lines[i][0]
            c_start = int(lst_refined_lines[i][1])
            c_end = int(lst_refined_lines[i][2])
            c_prob = float(lst_refined_lines[i][3])

            if (i != (len(lst_refined_lines) - 1)) and (lst_refined_lines[i + 1][0] == lst_refined_lines[i][0]):
                n_start = int(lst_refined_lines[i + 1][1])
                n_end = int(lst_refined_lines[i + 1][2])
                n_prob = float(lst_refined_lines[i + 1][3])

                dist = n_start - c_end

                if dist <= args.margin:
                    c_end = n_end
                    c_prob = max(c_prob, n_prob)

                    lst_refined_lines[i + 1][1] = str(c_start)
                    n_start = c_start
                    lst_refined_lines[i + 1][3] = str(c_prob)
                    n_prob = c_prob
                else:
                    pass
            else:
                n_start = c_end

            length = c_end - c_start

            if n_start != c_start or length > args.margin:
                if (c_prob == 0) or (c_prob == 1):
                    c_prob = int(c_prob)
                else:
                    pass

                new_refined.write(c_chrom + "\t" + str(c_start) + "\t" + str(c_end) + "\t" + str(c_prob) + "\n")
            else:
                pass
else:
    pass