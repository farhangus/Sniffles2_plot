# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC
"""
from arg_parser import argument_parser
from functions_variables import allele_frequency_chart_genrator
from functions_variables import samples_sv_numbers
def main():
    input_vcf_file,output_chrt=argument_parser()
#    allele_frequency_chart_genrator(input_vcf_file,output_chrt)
    samples_sv_numbers(input_vcf_file,output_chrt)
# Entry point of the program
if __name__ == "__main__":
    main()

