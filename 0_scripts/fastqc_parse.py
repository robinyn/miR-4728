import argparse
import os.path

def init_args():
    # The argument parser is created and the arguments the user can/must provide are defined
    parser = argparse.ArgumentParser(prog = "FastQC Parser", \
                                     description="Parses FastQC reports into a more R friendly form")
    parser.add_argument('-i', '--input', default="files_list.txt", nargs="?", \
                        help="Directory and name of the input file. A text (.txt) file with a single file on each line should be used for multiple files. (Default = files_list.txt)")
    parser.add_argument('-o', '--output', default="output.tsv", nargs="?", \
                        help="Directory and name of the output file. (Default = output.tsv)")

    # The arguments that the user provided are parsed and stored into the appropriate variables
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output

    if not os.path.isfile(input_file):
        print("Invalid input file. Aborting...")
        exit()

    if os.path.isfile(output_file):
        print("WARNING! The specified output file already exists. The script will overwrite the file.")

    # The parsed arguments are returned
    return input_file, output_file

def generate_summary(input_file):
    summary_list=[["sample_name", "total_read_count", "gc_content", "min_read_len", "max_read_len", "mean_read_len"]]

    try:
        with open(input_file, "r") as file_list:
            for file in file_list:
                report_dir = file.strip()

                total_read_count = 0
                gc_content = 0
                min_read_length = 0
                max_read_length = 0
                mean_read_length = 0

                length_module = False
                cur_len = 0
                len_count = 0
                len_sum = 0

                if report_dir.startswith("~"):
                    report_dir = os.path.expanduser(report_dir)

                with open(report_dir, "r") as report_file:
                    for line in report_file:
                        if line.startswith("Filename"):
                            sample_name = line.strip().split("\t")[1]

                        if line.startswith("Total Sequences"):
                            total_read_count = int(line.strip().split("\t")[1])
                            next

                        if line.startswith("Sequence length"):
                            min_read_length = int(line.strip().split("\t")[1].split("-")[0])
                            max_read_length = int(line.strip().split("\t")[1].split("-")[1])
                            next

                        if line.startswith("%GC"):
                            gc_content = int(line.strip().split("\t")[1])
                            next

                        if line.startswith("#length"):
                            length_module = True
                            next

                        if length_module and not line.startswith(">>END_MODULE"):
                            cur_len = int(line.strip().split("\t"[0]))
                            len_count = float(line.strip().split("\t"[1]))

                            len_sum += cur_len * len_count
                            next
                        elif length_module and line.startswith(">>END_MODULE"):
                            length_module = False
                            mean_read_length = len_sum / total_read_count
                            break
                summary_list.append([sample_name, str(total_read_count), str(gc_content), str(min_read_length), str(max_read_length), str(mean_read_length)])

    except Exception as e:
        print(e)

    return(summary_list)

def export_file(output_file, summary_list):
    try:
        with open(output_file, "w") as output_stream:
            for line in summary_list:
                print(line)
                output_stream.write("\t".join(line) + "\n")
    except Exception as e:
        print(e)



input_file, output_file = init_args()
summary_file = generate_summary(input_file)
export_file(output_file, summary_file)