# ===================================================================================================================================================================
# Title: novoalign_summaryParser.py
# Author: Euisuk Robin Han
# Description: A parser for the alignment summary from Novoalign.
# Date: 17/Apr/23
# ===================================================================================================================================================================

import argparse
import os.path

def init_args():
    # The argument parser is created and the arguments the user can/must provide are defined
    parser = argparse.ArgumentParser(prog = "HISAT2 Summary Parser", \
                                     description="Script to parse multiple alignment summary files from HISAT2 into one TSV file. Currently only works with paired reads.")
    parser.add_argument('-i', '--input', default="input.txt", nargs="?", \
                        help="Directory and name of the input file. A text (.txt) file with a single file on each line should be used for multiple files. (Default = input.sam.summary)")
    parser.add_argument('-o', '--output', default="parsed_output.tsv", nargs="?", \
                        help="Directory and name of the output file. (Default = parsed_ouput.tsv)")
    parser.add_argument('-e', '--extension', default=".sam.summary", nargs="?", \
                        help="Extensions/suffix to remove from the file names to extract the basename for the alignment samples. (Default = .sam.summary)")

    # The arguments that the user provided are parsed and stored into the appropriate variables
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    file_extension = args.extension

    if not os.path.isfile(input_file):
        print("Invalid input file. Aborting...")
        exit()

    if os.path.isfile(output_file):
        print("WARNING! The specified output file already exists. The script will overwrite the file.")

    # The parsed arguments are returned
    return input_file, output_file, file_extension

def parse_summary(input_file, output_file, file_extension):

    try:
        with open(input_file, "r") as inFile, open(output_file, "w") as outFile:
            outFile.write("sampleID\ttot_seq_count\tuniq_count\tmulti_count\tno_count\n")
            for file_name in inFile:
                sample_name = file_name.strip().split("/")[-1].removesuffix(file_extension)
                print(sample_name)
                with open(file_name.strip(), "r") as summary_file:
                    for line in summary_file:
                        if "Read Sequences" in line:
                            seq_count = line.split(":")[-1].strip()
                            print(seq_count)
                            continue
                        elif "Unique Alignment" in line:
                            uniq_map = line.split(":")[-1].strip().split(" ")[0].strip()
                            print(uniq_map)
                            continue
                        elif "Multi Mapped" in line:
                            multi_map = line.split(":")[-1].strip().split(" ")[0].strip()
                            print(multi_map)
                            continue
                        elif "No Mapping Found" in line:
                            no_map = line.split(":")[-1].strip().split(" ")[0].strip()
                            print(no_map)
                            break
                outFile.write("{}\t{}\t{}\t{}\t{}\n".format(sample_name, seq_count, uniq_map, multi_map, no_map))

    except Exception as e:
        print(e)

args = init_args()
parse_summary(args[0], args[1], args[2])