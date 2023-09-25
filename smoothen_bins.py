import argparse

parser = argparse.ArgumentParser(description = "Smoothens inputted bins using averages in rolling windows with user-configured sizes.")
parser.add_argument("bins", type = str, help = "Path to the text file that contains all the bins along with their scores. \
                    Assumed format for the file: CHROMOSOME START_COORDINATE   END_COORDINATE  SCORE/PROBABILITY")
parser.add_argument("output", type = str, help = "Path to the output file (text). The input format will be preserved.")
parser.add_argument("--window_size", type = int, default = 5, help = "Size of the rolling window - 5 by default. Has to be an ODD NUMBER.")
args = parser.parse_args()

if args.window_size % 2 == 0:
    print("Rolling window size must be odd!")
else:
    try:
        with open(args.bins, "r") as bins:
            bins_lines = bins.readlines()
    except FileNotFoundError:
        print("The bin file does not exist.")
        exit()

    bins_list = []

    for line in bins_lines:
        if line[0] != "#":
            elts = line.strip().split()
            ord_lst = [int(elts[1]), int(elts[2]), float(elts[3]), elts[0]]
            bins_list.append(ord_lst)
        else:
            pass

    ptr = 0

    while ptr < (len(bins_list) - (args.window_size // 2)):    
        if ptr == 0:
            check = 0
            for i in range(0, args.window_size):
                if bins_list[ptr + i][3] == bins_list[ptr][3]:
                    check += 1
                else: 
                    pass

            if (check == args.window_size):
                probs = [bins_list[j][2] for j in range(ptr, ptr + args.window_size)]
                local_sum = sum(probs)
                local_avg = local_sum / args.window_size
                bins_list[ptr + (args.window_size // 2)][2] = local_avg
            else:
                pass

            ptr += ((args.window_size // 2) + 1)
        else:
            check = 0
            for i in range(-(args.window_size // 2), (args.window_size // 2) + 1):
                if bins_list[ptr + i][3] == bins_list[ptr][3]:
                    check += 1
                else: 
                    pass

            if (check == args.window_size):
                probs = [bins_list[j][2] for j in range(ptr - (args.window_size // 2), ptr + (args.window_size // 2) + 1)]
                local_sum = sum(probs)
                local_avg = local_sum / args.window_size
                bins_list[ptr][2] = local_avg
            else:
                pass
            
            ptr += 1

    with open(args.output, "w") as output_file:
        output_file.write("# CHROMOSOME START_COORDINATE   END_COORDINATE  SCORE/PROBABILITY\n")
        
        for bin in bins_list:
            output_file.write(bin[3] + "\t" + str(bin[0]) + "\t" + str(bin[1]) + "\t" + str(bin[2]) + "\n")