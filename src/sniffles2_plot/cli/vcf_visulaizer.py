# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC - Farhang Jaryani
"""
import sys
sys.path.append('src/')  
from arg_parser import argument_parser
from sniffles2_plot.cli import single_visulaizer, multi_visulaizer
import os
def multi_vcf_flag(input_vcf_file):
    with open(input_vcf_file, "r") as f:
        for line in f:
            if line[0] != "#":
                if line.count('\t') > 9:
                    return 1
    return 0

def main():
    input_vcf_file,output_path=argument_parser()
    if (multi_vcf_flag(input_vcf_file)):
        multi_visulaizer(input_vcf_file, output_path)
    else:
        single_visulaizer(input_vcf_file, output_path)

if __name__ == "__main__":
    main()
