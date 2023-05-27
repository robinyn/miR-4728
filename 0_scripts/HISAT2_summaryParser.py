# ===================================================================================================================================================================
# Title: HISAT2_summaryParser.py
# Author: Euisuk Robin Han
# Description: A parser for the alignment summary from HISAT2. Currently only works with paired reads
# Date: 28/Mar/23
# ===================================================================================================================================================================

import argparse
import os.path

def init_args():
    # The argument parser is created and the arguments the user can/must provide are defined
    parser = argparse.ArgumentParser(prog = "HISAT2 Summary Parser", \
                                     description="Script to parse multiple alignment summary files from HISAT2 into one TSV file. Currently only works with paired reads.")
    parser.add_argument('-i', '--input', default="input.sam.summary", nargs="?", \
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
    summary_dict = {}
    num_samples = 0
    # Open summary files and parse results
    try:
        with open(input_file, "r") as inFile, open(output_file, "w") as outFile:
            if input_file.split(".")[-1] != "txt":
                basename = input_file.split("/")[-1].removesuffix(file_extension)

                print(basename)

                summary_dict["SampleID"]=basename
                for line in inFile:
                    if line.startswith("HISAT2"):
                        continue

                    line = line.strip()
                    category = line.split(":")[0].strip()
                    count = line.split(":")[1].strip().split(" ")[0]

                    if category == "Overall alignment rate":
                        count = count.removesuffix("%")

                    summary_dict[category]=count

                outFile.write("\t".join(summary_dict.keys()) + "\n")
                outFile.write("\t".join(summary_dict.values()))

            else:
                for index, file_name in enumerate(inFile):

                    if file_name == "":
                        continue

                    with open(file_name.strip(), "r") as summary_files:
                        basename = file_name.strip().split("/")[-1].removesuffix(file_extension)
                        print(basename)
                        summary_dict["SampleID"]=basename
                        for line in summary_files:
                            if line.startswith("HISAT2"):
                                continue

                            line = line.strip()
                            category = line.split(":")[0].strip()
                            count = line.split(":")[1].strip().split(" ")[0]

                            if category == "Overall alignment rate":
                                count = count.removesuffix("%")

                            summary_dict[category]=count

                    # Output to file
                    if index==0:
                        outFile.write("\t".join(summary_dict.keys()) + "\n")
                        outFile.write("\t".join(summary_dict.values()))
                    else:
                        outFile.write("\n" + "\t".join(summary_dict.values()))

    except Exception as e:
        print(e)

args = init_args()
parse_summary(args[0], args[1], args[2])
