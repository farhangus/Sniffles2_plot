# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC - Farhang Jaryani
"""
from sniffles2_plot.cli import generate_multi_vcf_charts, single_visulaizer

HEADER_SIGN = "#"


def _is_multi_vcf(input_vcf_file_path: str) -> bool:
    with open(input_vcf_file_path, "r") as f:
        for line in f:
            if not line.startswith(HEADER_SIGN):
                if line.count("\t") > 9:
                    return True
    return False


def generate_charts(input_vcf_file_path: str, output_directory_path: str) -> None:
    if _is_multi_vcf(input_vcf_file_path):
        generate_multi_vcf_charts(input_vcf_file_path, output_directory_path)
    else:
        single_visulaizer(input_vcf_file_path, output_directory_path)
