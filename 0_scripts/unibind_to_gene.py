# ===================================================================================================================================================================
# Title: unibind_to_gene.py
# Author: Euisuk Robin Han
# Description: A parser for mapping unibind data to genes
# Date: 02/May/23
# ===================================================================================================================================================================

import argparse
import os.path
import csv
from tqdm import tqdm

def init_args():
    # The argument parser is created and the arguments the user can/must provide are defined
    parser = argparse.ArgumentParser(prog = "Genome annotation parser", \
                                     description="A parser for genome annotation GTF files.")
    parser.add_argument('-g', '--gtf', default="input.gtf", nargs="?", \
                        help="Directory and name of the input GTF file.")
    parser.add_argument('-u', '--unibind', default="input.txt", nargs="?", \
                        help="Directory and name of the input GTF file.")
    parser.add_argument('-o', '--output', default="output.tsv", nargs="?", \
                        help="Directory and name of the output file.")
    parser.add_argument('-f', '--feature', default="gene", nargs="?", \
                        help="Genomic feature to parse")
    parser.add_argument('-a', '--attribute', default="gene_id", nargs="?", \
                        help="Annotation attribute to parse")

    # The arguments that the user provided are parsed and stored into the appropriate variables
    args = parser.parse_args()
    gtf_file = args.gtf
    unibind_file = args.unibind
    output_file = args.output
    feature = args.feature
    attribute = args.attribute

    if not os.path.isfile(gtf_file):
        print("ERROR! Invalid input GTF file. Aborting...")
        exit()

    if not os.path.isfile(unibind_file):
        print("ERROR! Invalid input UniBind file. Aborting...")
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
    return gtf_file, unibind_file, output_file, feature, attribute

def parse_gtf(input_file, parse_feature, parse_attr):
    parsed_anot = []
    chromosome_index = {}
    chromosome = ""
    feature = ""
    start_pos = ""
    end_pos = ""
    strand = ""
    attribute = ""

    try:
        num_lines = sum(1 for line in open(input_file, "r"))
        with open(input_file, "r") as gtf_file:
            for line in tqdm(gtf_file, total=num_lines):
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

                if "chr" not in chromosome:
                    continue

                attr_list = line[8].split(";")

                for attr in attr_list:
                    if attr.strip().startswith(parse_attr):
                        attribute = attr.strip().split(" ")[1].removeprefix('"').removesuffix('"')
                        continue

                parsed_anot.append([attribute, feature, chromosome, strand, start_pos, end_pos])

                if chromosome in chromosome_index.keys():
                    chromosome_index[chromosome].append(len(parsed_anot) - 1)
                else:
                    chromosome_index[chromosome] = [len(parsed_anot) - 1]


    except Exception as e:
        print("Error {}".format(e))

    return parsed_anot, chromosome_index

def parse_unibind(input_file):
    parsed_unibind = []
    chromosome = ""
    start_pos = ""
    end_pos = ""
    strand = ""
    cell_type = ""
    antibody = ""

    try:
        with open(input_file, "r") as unibind_file:
            num_lines = sum(1 for line in open(input_file, "r"))
            for line in tqdm(unibind_file, total=num_lines):
                if line.startswith("filename"):
                    continue

                line = line.strip().split("\t")

                cell_type = line[1]
                antibody = line[5]
                chromosome = line[11]
                start_pos = line[12]
                end_pos = line[13]
                strand = line[15]

                if "_" in chromosome:
                    continue

                parsed_unibind.append([cell_type, antibody, chromosome, strand, start_pos, end_pos])
    except Exception as e:
        print(e)

    return parsed_unibind

def coordinate_to_gene(annotation, chromosome_index, antibody, chromosome, strand, start_pos, end_pos):

    for index in chromosome_index[chromosome]:
        entry = annotation[index]

        if strand != entry[3]:
            continue

        if (int(start_pos) > int(entry[4])-1000) and (int(end_pos) < int(entry[4])+200):
            # print("{} is a target for {}".format(entry[0], antibody))
            # print("Promoter region: {} - {}".format(int(entry[4])-1000, int(entry[4])+200))
            # print("TFBS {} - {}".format(start_pos, end_pos))
            unibind_gene.append([antibody, entry[0]])

print("\nInitializing arguments...")
gtf_file, unibind_file, output_file, feature, attribute = init_args()

print("\nParsing gene annotation file...")
gencode, chromnosome_index = parse_gtf(gtf_file, feature, attribute)

print("\nParsing UniBind file...")
unibind = parse_unibind(unibind_file)

print("Converting coordinates to genes")
unibind_gene = []
for entry in tqdm(unibind, len(unibind_gene)):
    coordinate_to_gene(gencode, chromnosome_index, entry[1], entry[2], entry[3], entry[4], entry[5])

print("Writing output")

try:
    with open("unibind_genes.txt", "w") as output_stream:
        writer = csv.writer(output_stream, delimiter="\t")
        writer.writerows(unibind_gene)
except Exception as e:
    print(e)
