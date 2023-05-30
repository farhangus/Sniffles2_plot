# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""
import subprocess
from arg_parser import argument_parser
from functions_variables_multi import *

def main():
    input_vcf_file,output_chrt=argument_parser()
    Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_chrt)
    Genome_chart_data_generator.allele_frequency_chart_generator()
    Genome_chart_data_generator.samples_sv_numbers()
    cmd= f"rm -f  {output_chrt}*.txt"
    subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    main()

