# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 10:48:59 2023

@author: HGSC
"""

import sys
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description="sample1111")
    parser.add_argument("-i", "--input", help="Specify the path to the VCF files directory (accepts both single and multi samples VCF files)",required=True)
    parser.add_argument("-o", "--output", help="Specify the path to the VCF outp file")
    args = parser.parse_args()
    input_file_path,output_file_path=args.input,args.output
    if os.path.isfile(input_file_path):
            file_name = os.path.splitext(input_file_path)[0]
            print(file_name)
            directory = os.path.dirname(os.path.abspath(input_file_path))
            os.makedirs(directory, exist_ok=True)
            print("Created directory:", directory)
            os.system(f"python3 vcf_visulaizer.py -i {input_file_path} -o {output_file_path}/")
    else:
        for file_entry in os.scandir(input_file_path):
            file_name=file_entry.name
            if file_name.lower().endswith(".vcf"):
                file_path = file_entry.path
                if os.path.isfile(file_path):
                    file_name = os.path.splitext(file_name)[0]
                    directory_path = os.path.join(input_file_path, file_name)
                    os.makedirs(directory_path, exist_ok=True)
                    print("Created directory:", directory_path)
                    os.system(f"python3 vcf_visulaizer.py -i {file_path} -o {directory_path}/")
if __name__ == "__main__":
    main()
