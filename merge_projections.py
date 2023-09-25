import argparse
import sys

def any(lst, string):
    for item in lst:
        if item in string:
            return True
        else:
            pass
    return False

def find_shortest(entry_list, completed_entries):
    try:
        min_length = min([e[2] for e in entry_list if e[2] != 0])
    except ValueError:
        min_length = 0

    for i, entry in enumerate(entry_list):
        if (min_length != 0) and (entry[2] == min_length) and (entry not in completed_entries):
            return i, entry
        else:
            pass
    
    return sys.maxsize, None

def calc_length(entry):
    return (int(entry[0][1]) - int(entry[0][0]))

parser = argparse.ArgumentParser(description = "Merges genomic region projections using the user-provided margin as basis.")
parser.add_argument("input", type = str, help = "Path to the list of projected genomic regions. Please use a header to describe the columns. \
                    Assumed format for the file: CHROMOSOME  START_COORDINATE    END_COORDINATE  ACTIVE_TISSUES  TISSUE_1_SCORE  [TISSUE_2_SCORE] ... \
                    The 4th column MUST be named 'ACTIVE_TISSUES,' which is anyway the case in the output of project_regions.py. \
                    Columns containing the tissue scores 1) MUST directly follow 'ACTIVE_TISSUES', and 2) contain the tissue name within their name (e.g. HL_SCORE, where HL refers to a tissue), \
                    although the order of tissues does not matter.")
parser.add_argument("--output", type = str, default = "stdout", help = "Path to the output file (text). \
                    If specified as not a proper name but 'stdout', which is the default option, the output will simply be printed and can be piped.")
parser.add_argument("--margin", type = int, default = 100, help = "The minimum length allowed for a region. Regions smaller than or equal to the \
                    margin in size will either be split between or merged directly with the preceding and/or succeeding region depending on the specific case, \
                    taking scores and active tissues into consideration. 100 by default.")
parser.add_argument("--sep", type = str, default = ",", help = "The separator for tissue names under 'ACTIVE_TISSUES'. ',' by default.")
parser.add_argument("--len", action = 'store_true', help = "The final length of each region will be included as an extra column in the output if this option is provided.")
args = parser.parse_args()

region_dict = []
header = []
dict_template = []
all_tissues = set()

try:
    with open(args.input, "r") as input:
        input_lines = input.readlines()
except FileNotFoundError: 
    print("The input file doesn't exist.")
    exit()

for i, line in enumerate(input_lines):
    if line[0] != "#":
        elts = line.strip().split()
        
        coords = [int(elts[1]), int(elts[2]), elts[0]]
        dict_construct = {}
        length = int(elts[2]) - int(elts[1])

        region_dict.append([coords, {}, length])

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

        for tissues in [region_dict[i - 1][1]["ACTIVE_TISSUES"]]:
            for tissue in tissues.strip().split(args.sep):
                all_tissues.add(tissue)
    else:
        header = line.strip().split()
        
        for i, col in enumerate(header):
            if i != 0 and i != 1 and i != 2 and i != 3:
                dict_template.append(col)
            else:
                pass

blocks = []
block = []

for i in range(0, len(region_dict)):
    if i != len(region_dict) - 1 and region_dict[i][0][1] == region_dict[i + 1][0][0]:
        block.append(region_dict[i])
    else:
        block.append(region_dict[i])
        blocks.append(block)
        block = []

tissue_count = len(all_tissues)

for block in blocks:
    completed_entries = []
    i, shortest_entry = find_shortest(block, completed_entries)

    while shortest_entry != None:
        if i != 0:
            prev = i - 1
        else:
            prev = 0

        if i != len(block) - 1:
            next = i + 1
        else:
            next = 0
        
        while prev != 0 and block[prev][1] == {}:
            prev -= 1

        while next < len(block) and block[next][1] == {}:
            next += 1

        if shortest_entry[2] <= args.margin:
            if i != 0:
                if block[prev][1] != {}:
                    prev_tissues = block[prev][1]["ACTIVE_TISSUES"].strip().split(args.sep)
                else:
                    prev_tissues = []
            else:
                prev_tissues = []
            
            current_tissues = shortest_entry[1]["ACTIVE_TISSUES"].strip().split(args.sep)

            if i != len(block) - 1:
                if block[next][1] != {}:
                    next_tissues = block[next][1]["ACTIVE_TISSUES"].strip().split(args.sep)
                else:
                    next_tissues = []
            else:
                next_tissues = []
            
            if all(tissue in prev_tissues for tissue in current_tissues) == True and all(tissue in next_tissues for tissue in current_tissues) == True and block[prev][0][2] == shortest_entry[0][2] == block[next][0][2] and block[prev][0][1] == shortest_entry[0][0] and block[next][0][0] == shortest_entry[0][1]: # SPLITTING
                if all(tissue in prev_tissues for tissue in next_tissues) == True:
                    block[next][0][0] = block[prev][0][0]

                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        block[next][1][key] = max(block[prev][1][key], shortest_entry[1][key], block[next][1][key])

                    block[i] = [[], {}, 0]
                    block[prev] = [[], {}, 0]

                    block[next][2] = calc_length(block[next])
                else:
                    region_length = shortest_entry[2]
                    block[next][0][0] -= (region_length / 2)
                    block[prev][0][1] += (region_length / 2)

                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        block[next][1][key] = max(block[next][1][key], shortest_entry[1][key])
                        block[prev][1][key] = max(block[prev][1][key], shortest_entry[1][key])

                    block[i] = [[], {}, 0]

                    block[next][2] = calc_length(block[next])
                    block[prev][2] = calc_length(block[prev])
            elif all(tissue in prev_tissues for tissue in current_tissues) == True and block[prev][0][2] == shortest_entry[0][2] and block[prev][0][1] == shortest_entry[0][0]: # DIRECT ASSIGNMENT -> PREVIOUS
                block[i][0][0] = block[prev][0][0]
                
                for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                    block[i][1][key] = max(block[prev][1][key], shortest_entry[1][key])

                block[i][1]["ACTIVE_TISSUES"] = block[prev][1]["ACTIVE_TISSUES"]

                block[prev] = [[], {}, 0]

                block[i][2] = calc_length(block[i])
            elif all(tissue in next_tissues for tissue in current_tissues) == True and block[next][0][2] == shortest_entry[0][2] and block[next][0][0] == shortest_entry[0][1]: # DIRECT ASSIGNMENT -> NEXT
                block[next][0][0] = shortest_entry[0][0]
                
                for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                    block[next][1][key] = max(block[next][1][key], shortest_entry[1][key])
                
                block[i] = [[], {}, 0]

                block[next][2] = calc_length(block[next])
            elif all(tissue in prev_tissues for tissue in current_tissues) == False and all(tissue in next_tissues for tissue in current_tissues) == False:
                if (any(current_tissues, prev_tissues) == True and any(current_tissues, next_tissues) == True) and (block[prev][0][2] == shortest_entry[0][2] == block[next][0][2] and block[prev][0][1] == shortest_entry[0][0] and block[next][0][0] == shortest_entry[0][1]):
                    score_change_next = 0
                    score_change_prev = 0

                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        if any(shortest_entry[1]["ACTIVE_TISSUES"].strip().split(args.sep), key) == True:
                            score_change = block[next][1][key] - shortest_entry[1][key]

                            if score_change > score_change_next:
                                score_change_next = score_change
                            else:
                                pass
                        else:
                            pass
                    
                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        if any(shortest_entry[1]["ACTIVE_TISSUES"].strip().split(args.sep), key) == True:
                            score_change = block[prev][1][key] - shortest_entry[1][key]

                            if score_change > score_change_prev:
                                score_change_prev = score_change
                            else:
                                pass
                        else:
                            pass
                    
                    if score_change_prev > score_change_next: # SCORE-BASED ASSIGNMENT -> PREVIOUS
                        block[i][0][0] = block[prev][0][0]
                    
                        for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                            block[i][1][key] = max(block[prev][1][key], shortest_entry[1][key])
                        
                        if block[prev][2] > shortest_entry[2]:
                            block[i][1]["ACTIVE_TISSUES"] = block[prev][1]["ACTIVE_TISSUES"]
                        else:
                            pass

                        block[prev] = [[], {}, 0]

                        block[i][2] = calc_length(block[i])
                    elif score_change_next > score_change_prev: # SCORE-BASED ASSIGNMENT -> NEXT
                        block[next][0][0] = shortest_entry[0][0]
                    
                        for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                            block[next][1][key] = max(block[next][1][key], shortest_entry[1][key])
                        
                        if shortest_entry[2] > block[next][2]:
                            block[next][1]["ACTIVE_TISSUES"] = shortest_entry[1]["ACTIVE_TISSUES"]
                        else:
                            pass

                        block[i] = [[], {}, 0]

                        block[next][2] = calc_length(block[next])
                    elif score_change_next == score_change_prev: # SPLITTING
                        if all(tissue in prev_tissues for tissue in next_tissues) == True:
                            block[next][0][0] = block[prev][0][0]

                            for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                                block[next][1][key] = max(block[prev][1][key], shortest_entry[1][key], block[next][1][key])

                            max_length = max(block[prev][2], shortest_entry[2], block[next][2])

                            if block[next][2] == max_length:
                                pass
                            elif block[prev][2] == max_length:
                                block[next][1]["ACTIVE_TISSUES"] = prev_tissues
                            elif shortest_entry[2] == max_length:
                                block[next][1]["ACTIVE_TISSUES"] = current_tissues
                            else:
                                pass

                            block[i] = [[], {}, 0]
                            block[prev] = [[], {}, 0]

                            block[next][2] = calc_length(block[next])
                        else:
                            region_length = shortest_entry[2]
                            block[next][0][0] -= (region_length / 2)
                            block[prev][0][1] += (region_length / 2)

                            for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                                block[next][1][key] = max(block[next][1][key], shortest_entry[1][key])
                                block[prev][1][key] = max(block[prev][1][key], shortest_entry[1][key])
                            
                            max_length = max(block[prev][2], (shortest_entry[2] / 2), block[next][2])

                            if block[next][2] == max_length:
                                pass
                            elif block[prev][2] == max_length:
                                pass
                            elif (shortest_entry[2] / 2) == max_length:
                                block[next][1]["ACTIVE_TISSUES"] = shortest_entry[1]["ACTIVE_TISSUES"]
                                block[prev][1]["ACTIVE_TISSUES"] = shortest_entry[1]["ACTIVE_TISSUES"]
                            else:
                                pass

                            block[i] = [[], {}, 0]

                            block[next][2] = calc_length(block[next])
                            block[prev][2] = calc_length(block[prev])
                    else:
                        pass
                elif any(current_tissues, prev_tissues) == True and block[prev][0][2] == shortest_entry[0][2] and block[prev][0][1] == shortest_entry[0][0]: # DIRECT ASSIGNMENT -> PREVIOUS
                    block[i][0][0] = block[prev][0][0]
                    
                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        block[i][1][key] = max(block[prev][1][key], shortest_entry[1][key])
                    
                    if block[prev][2] > shortest_entry[2]:
                        block[i][1]["ACTIVE_TISSUES"] = block[prev][1]["ACTIVE_TISSUES"]
                    else:
                        pass

                    block[prev] = [[], {}, 0]

                    block[i][2] = calc_length(block[i])
                elif any(current_tissues, next_tissues) == True and block[next][0][2] == shortest_entry[0][2] and block[next][0][0] == shortest_entry[0][1]: # DIRECT ASSIGNMENT -> NEXT
                    block[next][0][0] = shortest_entry[0][0]
                    
                    for key in list(shortest_entry[1].keys())[1:1 + tissue_count + 1]:
                        block[next][1][key] = max(block[next][1][key], shortest_entry[1][key])
                    
                    if shortest_entry[2] > block[next][2]:
                        block[next][1]["ACTIVE_TISSUES"] = shortest_entry[1]["ACTIVE_TISSUES"]
                    else:
                        pass

                    block[i] = [[], {}, 0]

                    block[next][2] = calc_length(block[next])
                else:
                    block[i] = [[], {}, 0]
            else:
                block[i] = [[], {}, 0]
        else:
            pass

        if ((block[i] == [[], {}, 0]) or (block[prev] == [[], {}, 0]) or (block[next] == [[], {}, 0])) and ([[], {}, 0] not in completed_entries):
            completed_entries.append([[], {}, 0])
        else:
            pass
        
        if shortest_entry[2] > args.margin:
            completed_entries.append(shortest_entry)
        else:
            pass

        i, shortest_entry = find_shortest(block, completed_entries)

if args.len == True:
    header.append("LENGTH")
else:
    pass

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

    for block in blocks:
        for entry in block:
            if entry[1] != {}:
                print(entry[0][2] + "\t" + str(int(entry[0][0])) + "\t" + str(int(entry[0][1])), end = "")

                for i, key in enumerate(list(entry[1].keys())):
                    if (entry[1][key] == 0) or (entry[1][key] == 1):
                        entry[1][key] = int(entry[1][key])
                    else:
                        pass

                    if i == (len(list(entry[1].keys())) - 1):
                        if args.len == True:
                            print("\t" + str(entry[1][key]), end = "")
                            print("\t" + str(entry[2]))
                        else:
                            print("\t" + str(entry[1][key]))
                    else:
                        print("\t" + str(entry[1][key]), end = "")
            else:
                pass
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

        for block in blocks:
            for entry in block:
                if entry[1] != {}:
                    output.write(entry[0][2] + "\t" + str(int(entry[0][0])) + "\t" + str(int(entry[0][1])))

                    for i, key in enumerate(list(entry[1].keys())):
                        if (entry[1][key] == 0) or (entry[1][key] == 1):
                            entry[1][key] = int(entry[1][key])
                        else:
                            pass
                        
                        if i == (len(list(entry[1].keys())) - 1):
                            if args.len == True:
                                output.write("\t" + str(entry[1][key]))
                                output.write("\t" + str(entry[2]) + "\n")
                            else:
                                output.write("\t" + str(entry[1][key]) + "\n")
                        else:
                            output.write("\t" + str(entry[1][key]))
                else:
                    pass