#! /usr/bin/bash

files=()
bin_files=()
naming_symbol="_"
run_script=false
keep_temp=false
margin="0"
output=""
final_output="stdout"
project_path="$(dirname "$0")/project_regions.py"
merge_path="$(dirname "$0")/merge_projections.py"
recompute_path="$(dirname "$0")/recompute_scores.py"
keep_path=""

while [[ $# -gt 0 ]]; do
    key="$1"
    
    case $key in
        -f|--files)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -f!"
                exit
            else
                while [[ $# -gt 0 && "$2" != -* ]]; do
                    files+=("$2")
                    shift
                done
                shift
            fi
        ;;
        -bf|--bin_files)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -o!"
                exit
            else
                while [[ $# -gt 0 && "$2" != -* ]]; do
                    bin_files+=("$2")
                    shift
                done
                shift
            fi
        ;;
        -ns|--naming_symbol)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -ns!"
                exit
            else
                naming_symbol="$2"
                shift
                shift
            fi
        ;;
        -o|--output)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -o!"
                exit
            else
                output="$2"
                shift
                shift
            fi
        ;;
        -r|--run)
            run_script=true
            shift
        ;;
        -fo|--final_output)
            if [ "$run_script" = false ]; then
                echo "You have to run the script first to set a final output! Add -r to your command."
                exit
            elif [ -z "$2" ]; then
                echo "You have to provide an argument for -fo!"
                exit
            else
                final_output="$2"
                shift
                shift
            fi
        ;;
        -m|--margin)
            if [ "$run_script" = false ]; then
                echo "You have to run the script first to set a margin! Add -r to your command."
                exit
            elif [ -z "$2" ]; then
                echo "You have to provide an argument for -m!"
                exit
            else
                margin="$2"
                shift
                shift
            fi
        ;;
        -projscript|--projection_script)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -projscript!"
                exit
            else
                project_path="$2"
                shift
                shift
            fi
        ;;
        -recompscript|--recomputation_script)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -recompscript!"
                exit
            else
                recompute_path="$2"
                shift
                shift
            fi
        ;;
        -mergscript|--merging_script)
            if [ -z "$2" ]; then
                echo "You have to provide an argument for -mergscript!"
                exit
            else
                merge_path="$2"
                shift
                shift
            fi
        ;;
        -k|--keep)
            if [ ! -z "$2" ]; then
                keep_path="$2"
                shift
            fi
            keep_temp=true
            shift
        ;;
        -h|--help)
            echo "Description: Wrapper script for project_regions.py and merge_projections.py."
            echo 
            echo "Usage: $0 [-f|--files] [-o|--output] [-ns|--naming_symbol] [-projscript|--projection_script] [-recompscript|--recomputation_script] [-bf|--bin_files] [-mergscript|---merging_script] [-r|--run] [-fo|--final_output] [-m|--margin] [-k|--keep]"
            echo 
            echo "REQUIRED: -f, -o"
            echo "OPTIONAL: -ns, -r, -fo, -m, -projscript, ((-recompscript, -bf), -mergscript), -k"
            echo
            echo "-f (FILES): Takes the files that contain the genomic regions to be inputted to the projection script. Regions in each file is labelled with the part of the file name that precedes the first naming symbol (e.g. with -ns _ (default naming symbol), regions in FL_CS17.merge.SMOOTHENED.singleEnh.bedGraph would be marked with FL)."
            echo "-o (OUTPUT): Takes the path to the output file, which will be a sorted concatenation of the provided genomic regions with their labels. Before concatenation and sorting, overlapping regions in each tissue are merged with each other."
            echo "-ns (NAMING_SYMBOL): The symbol used to parse tissue names. Text before the first occurrence of this symbol is recognised as the tissue name. '_' by default."
            echo "-r (RUN): Runs the Python script using the output file."
            echo "-fo (FINAL_OUTPUT): Takes the path to the final output file, that is, the output of the projection script. Works only when the -r option is provided (when the script is run). The final output is directly printed if not specified (allows for piping)."
            echo "-m (MARGIN): Takes the margin to be provided to the merging script. It refers to the minimum size allowed for a region. If a region is shorter than or equal to the margin in size, it will be merged with one of the surrounding regions depending on the specific case. 0 (turned off) by default."
            echo "-projscript (PROJECTION_SCRIPT): Takes the path to project_regions.py. If not specified, it is assumed that the Python script is found in the same directory as the wrapper script with the same name - project_regions.py (recommended)."
            echo "-recompscript (RECOMPUTATION_SCRIPT): Takes the path to recompute_scores.py. If not specified, it is assumed that the Python script is found in the same directory as the wrapper script with the name 'recompute_scores.py' (recommended). IT REQUIRES SPECIFYING THE PATHS TO THE BIN FILES THAT WILL BE USED FOR SCORE RECOMPUTATION WITH -bf."
            echo "-bf (BIN_FILES): Takes the paths to the bin files that will be used for running -recompscript."
            echo "-mergscript (MERGING_SCRIPT): Takes the path to merge_projections.py. If not specified, it is assumed that the Python script is found in the same directory as the wrapper script with the name 'merge_projections.py' (recommended)."
            echo "-k (KEEP): If provided, temporary files that are created while generating the main output are kept instead of removed (e.g. sorted and merged versions of the files). If no argument is given, temporary files are kept in the working directory; if a path is provided, they are moved there."
            exit
        ;;
        *)
            echo "Invalid option! Please check $0 -h|--help."
            exit
        ;;
    esac
done

if [ ${#files[@]} -eq 0 ]; then
    echo "-f is NOT an optional argument!"
    exit
fi

if [ "$output" == "" ]; then
    echo "-o is NOT an optional argument!"
    exit
fi

if [ ${#bin_files[@]} -eq 0 ] && [ "$margin" -ne 0 ]; then
    echo "You have to specify the relevant bin files for score recomputation if you want to merge projections! Specify bin files with -bf FILE_1 FILE_2 ... FILE_N."
    exit
fi

function sort_regions {
    if [ -z "$1" ]; then
        for file in "${files[@]}"; do
            sort -b -k1,1 -k2n,2n "$file" > "$file.sorted"
        done
    elif [ "$1" = "labelled_concatenated" ]; then
        sort -b -k1,1 -k2n,2n -k3n,3n "$output.unordered" > "$output"
        sed -i '1 s/^/# CHROMOSOME START_COORDINATE    END_COORDINATE  TISSUE\n/' "$output"

        if [ "$keep_temp" = false ] ; then
            rm "$output.unordered"
        else
            if [ "$keep_path" != "" ]; then
                if [ ! -d "$keep_path" ]; then
                    mkdir -p "$keep_path"
                fi
                mv "$output.unordered" "$keep_path"
            fi
        fi
    fi
}

function merge_regions {
    for file in "${files[@]}"; do
        bedtools merge -i "$file.sorted" -c 4 -o max > "$file.merged"
    done

    if [ "$keep_temp" = false ] ; then
        for file in "${files[@]}"; do
            rm "$file.sorted"
        done
    else
        if [ "$keep_path" != "" ]; then
            if [ ! -d "$keep_path" ]; then
                mkdir -p "$keep_path"
            fi
            for file in "${files[@]}"; do
                mv "$file.sorted" "$keep_path"
            done
        fi
    fi
}

function label_regions {
    for file in "${files[@]}"; do
        awk -F "\t" -v sep="$naming_symbol" 'BEGIN {OFS="\t"; filename=ARGV[1]; sub(".*/", "", filename); split(filename, a, sep); tissue=a[1]} {print $0, tissue;}' "$file.merged" > "$file.labelled"
    done

    if [ "$keep_temp" = false ] ; then
        for file in "${files[@]}"; do
            rm "$file.merged"
        done
    else
        if [ "$keep_path" != "" ]; then
            if [ ! -d "$keep_path" ]; then
                mkdir -p "$keep_path"
            fi
            for file in "${files[@]}"; do
                mv "$file.merged" "$keep_path"
            done
        fi
    fi
}

function concatenate_files {
    for file in "${files[@]}"; do
        cat "$file.labelled" >> "$output.unordered"
    done

    if [ "$keep_temp" = false ] ; then
        for file in "${files[@]}"; do
            rm "$file.labelled"
        done
    else
        if [ "$keep_path" != "" ]; then
            if [ ! -d "$keep_path" ]; then
                mkdir -p "$keep_path"
            fi
            for file in "${files[@]}"; do
                mv "$file.labelled" "$keep_path"
            done
        fi
    fi

    sort_regions labelled_concatenated
}

function run {
    if [ "$margin" -ne 0 ]; then
        python3 "$project_path" "$output" --output "$final_output.temp"
        python3 "$recompute_path" ${bin_files[*]} "$final_output.temp" --output "$final_output.temp_computed" --ns "$naming_symbol"
        python3 "$merge_path" "$final_output.temp_computed" --output "$final_output" --margin "$margin"

        if [ "$keep_temp" = false ] ; then
            rm "$final_output.temp"
            rm "$final_output.temp_computed"
        else
            if [ "$keep_path" != "" ]; then
                if [ ! -d "$keep_path" ]; then
                    mkdir -p "$keep_path"
                fi
                mv "$final_output.temp" "$keep_path"
                mv "$final_output.temp_computed" "$keep_path"
            fi
        fi
    else
        python3 "$project_path" "$output" --output "$final_output"
    fi
}

sort_regions
merge_regions
label_regions
concatenate_files

if [ "$run_script" = true ]; then
    run
fi
