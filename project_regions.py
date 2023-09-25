import argparse

def text_to_dict(input_name):
    region_dict = {} # FORMAT: (idx, chr): [[TISSUE_1, ..., TISSUE_N], [STATE_1, ..., STATE_N (+/-)]]

    try:
        with open(input_name, "r") as regions:
            regions_lines = regions.readlines()
    except FileNotFoundError:
        print("The input file doesn't exist.")
        exit()

    for line in regions_lines:
        if line[0] != "#":
            elts = line.strip().split()

            start_key = (int(elts[1]), elts[0])
            end_key = (int(elts[2]), elts[0])

            if start_key not in region_dict.keys():
                region_dict[start_key] = [[elts[4]], ["+"]]
            else:
                region_dict[start_key][0].append(elts[4])
                region_dict[start_key][1].append("+")

            if end_key not in region_dict.keys():
                region_dict[end_key] = [[elts[4]], ["-"]]
            else:
                region_dict[end_key][0].append(elts[4])
                region_dict[end_key][1].append("-")
        else:
            pass

    return region_dict

def dict_to_output(region_dict, output_name):
    sorted_keys = sorted(region_dict.keys(), key = lambda k: (k[1], k[0]))
    
    if output_name == "stdout":
        active_tissues = []
        start_idx = 0
        start_chrom = ""
        
        print("# CHROMOSOME    START_COORDINATE    END_COORDINATE  ACTIVE_TISSUES")

        for key in sorted_keys:
            for i in range(len(region_dict[key][0])):
                if (start_chrom != "" or start_chrom != key[1]) and (start_idx != 0 and start_idx != key[0]):
                    tissues = ",".join(sorted(list(set(active_tissues))))
                    print(start_chrom + "\t" + str(start_idx) + "\t" + str(key[0]) + "\t" + tissues)
                else:
                    pass
                
                if region_dict[key][1][i] == "+":
                    active_tissues.append(region_dict[key][0][i])
                    start_idx = key[0]
                    start_chrom = key[1]
                elif region_dict[key][1][i] == "-":
                    active_tissues.remove(region_dict[key][0][i])
                    
                    if len(active_tissues) == 0:
                        start_idx = 0
                        start_chrom = ""
                    else:
                        start_idx = key[0]
                        start_chrom = key[1]
                else:
                    pass
    else:
        with open(output_name, "w") as output:
            active_tissues = []
            start_idx = 0
            start_chrom = ""

            output.write("# CHROMOSOME    START_COORDINATE    END_COORDINATE  ACTIVE_TISSUES\n")

            for key in sorted_keys:
                for i in range(len(region_dict[key][0])):
                    if (start_chrom != "" or start_chrom != key[1]) and (start_idx != 0 and start_idx != key[0]):
                        tissues = ",".join(sorted(list(set(active_tissues))))
                        output.write(start_chrom + "\t" + str(start_idx) + "\t" + str(key[0]) + "\t" + tissues + "\n")
                    else:
                        pass
                    
                    if region_dict[key][1][i] == "+":
                        active_tissues.append(region_dict[key][0][i])
                        start_idx = key[0]
                        start_chrom = key[1]
                    elif region_dict[key][1][i] == "-":
                        active_tissues.remove(region_dict[key][0][i])
                        
                        if len(active_tissues) == 0:
                            start_idx = 0
                            start_chrom = ""
                        else:
                            start_idx = key[0]
                            start_chrom = key[1]
                    else:
                        pass

parser = argparse.ArgumentParser(description = "Projects given genomic regions over each other to reveal overlaps and non-overlaps.")
parser.add_argument("input", type = str, help = "Path to the MERGED, CONCATENATED and SORTED list of genomic regions to be projected on each other WITH THE TISSUE LABEL FOR EACH REGION. Refer to the wrapper script (project_regions.sh) if you are not sure about how to prepare the input file.")
parser.add_argument("--output", type = str, default = "stdout", help = "Path to the output file (text). If specified as not a proper name but 'stdout', which is the default option, the output will simply be printed and can be piped.")
args = parser.parse_args()

region_dict = text_to_dict(args.input)
dict_to_output(region_dict, args.output)