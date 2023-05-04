# ===================================================================================================================================================================
# Title: gtf_parser.py
# Author: Euisuk Robin Han
# Description: A parser for genome annotation GTF files to map features to genomic regions
# Date: 02/May/23
# ===================================================================================================================================================================

import argparse
import os.path

def init_args():
    # The argument parser is created and the arguments the user can/must provide are defined
    parser = argparse.ArgumentParser(prog = "Genome annotation parser", \
                                     description="A parser for genome annotation GTF files.")
    parser.add_argument('-i', '--input', default="input.gtf", nargs="?", \
                        help="Directory and name of the input file.")
    parser.add_argument('-o', '--output', default="output.tsv", nargs="?", \
                        help="Directory and name of the output file.")
    parser.add_argument('-f', '--feature', default="gene", nargs="?", \
                        help="Genomic feature to parse")
    parser.add_argument('-a', '--attribute', default="gene_id", nargs="?", \
                        help="Annotation attribute to parse")

    # The arguments that the user provided are parsed and stored into the appropriate variables
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output
    feature = args.feature
    attribute = args.attribute

    if not os.path.isfile(input_file):
        print("ERROR! Invalid input file. Aborting...")
        exit()

    if os.path.isfile(output_file):
        print("WARNING! The specified output file already exists. The script will overwrite the file.")

    if (feature=="gene") and (attribute=="transcript_id" or attribute=="exon_id"):
        print("ERROR! Cannot parse attribute '{}' for feature '{}'. Aborting...".format(attribute, feature))
        exit()

    if (feature=="transcript") and (attribute=="exon_id"):
        print("ERROR! Cannot parse attribute '{}' for feature '{}'. Aborting...".format(attribute, feature))
        exit()

    # The parsed arguments are returned
    return input_file, output_file, feature, attribute

def parse_gtf(input_file, output_file, parse_feature, parse_attr):
    parsed_anot = []
    chromosome = ""
    feature = ""
    start_pos = ""
    end_pos = ""
    strand = ""
    attribute = ""

    try:
        with open(input_file, "r") as gtf_file:
            for line in gtf_file:
                if line.startswith("##"):
                    continue
                line = line.split("\t")
                chromosome = line[0]
                feature = line[2]
                start_pos = line[3]
                end_pos = line[4]
                strand = line[6]

                if feature != parse_feature:
                    continue

                attr_list = line[8].split(";")

                for attr in attr_list:
                    if attr.strip().startswith(parse_attr):
                        attribute = attr.strip().split(" ")[1].removeprefix('"').removesuffix('"')
                        print(attribute)
                        continue

                parsed_anot.append([attribute, feature, chromosome, strand, start_pos, end_pos])

    except Exception as e:
        print(e)

    return parsed_anot

input_file, output_file, feature, attribute = init_args()
gencode = parse_gtf(input_file, output_file, feature, attribute)
gencode_dict = {}

for entry in gencode:
    gencode_dict[entry[0]] = entry[1:]

control = True

while control:
    gene = input("Input gene ID: ")
    if gene == "exit":
        control=False
        continue
    print(gencode_dict[gene])

