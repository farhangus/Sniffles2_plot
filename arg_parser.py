#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from argparse import RawTextHelpFormatter

def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__, prog='vcf_visulizer', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-i", "--input", dest="VCF_file",
                        help="input path to VCF file",required=True)
    parser.add_argument("-o", "--output", dest="output_graph",
                        help="output to the graph")



    args = parser.parse_args()
    
    input_file = args.VCF_file
    output_file = args.output_graph


    return input_file,output_file

