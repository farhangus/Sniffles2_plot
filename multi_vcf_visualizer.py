# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""
import subprocess
from arg_parser import argument_parser
from functions_variables_multi import *

def multi_visulaizer(input_vcf_file,output_chrt):
    #input_vcf_file,output_chrt=argument_parser()

    Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_chrt)
    Genome_chart_data_generator.allele_frequency_chart_generator()
    Genome_chart_data_generator.samples_sv_numbers()


if __name__ == "__main__":
    main()

