# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 10:48:59 2023

@author: HGSC
"""

import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Snfiffle2 plot generator")
    parser.add_argument(
        "-i",
        "--input",
        help="Specify the path to the VCF files directory (accepts both single and multi samples VCF files)",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="Specify the path to the VCF output file"
    )
    args = parser.parse_args()
    input_file_path, output_file_path = args.input, args.output
    if os.path.isfile(input_file_path):
        if os.path.exists(output_file_path) and not os.path.isdir(output_file_path):
            raise IOError("the out path is not a directory")
        if not os.path.exists(output_file_path):
            os.mkdir(output_file_path)
        file_name = os.path.splitext(input_file_path)[0]
        print(file_name)
        os.system(
            f"python3 src/sniffles2_plot/cli/vcf_visulaizer.py -i {input_file_path} -o {output_file_path}"
        )
    else:
        for file_entry in os.scandir(input_file_path):
            file_name = file_entry.name
            if file_name.lower().endswith(".vcf"):
                file_path = file_entry.path
                file_name = os.path.splitext(file_name)[0]
                directory_path = os.path.join(input_file_path, file_name)
                os.makedirs(directory_path, exist_ok=True)
                print("Created directory:", directory_path)
                os.system(
                    f"python3 src/sniffles2_plot/cli/vcf_visulaizer.py -i {file_path} -o {directory_path}"
                )


if __name__ == "__main__":
    main()
