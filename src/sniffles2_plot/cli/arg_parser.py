#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse


def argument_parser():
    """function to get the input and output file paths"""
    parser = argparse.ArgumentParser(
        description="A tool for visualizing VCF files.",
        prog="vcf_visualizer",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", dest="vcf_file", help="Input path to VCF file", required=True
    )
    parser.add_argument(
        "-o", "--output", dest="output_directory", help="Output path to the graph"
    )

    args = parser.parse_args()
    input_file = args.vcf_file
    output_directory = args.output_directory

    return input_file, output_directory
