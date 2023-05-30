# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC - Farhang Jaryani
"""
from arg_parser import argument_parser
from functions_variables import *
import os
sv_ranges = {
    "DEL": None,
    "INS": None,
    "DUP": None,
    "INV": None,
    "BND": None,
}

input_vcf_file,output_path=argument_parser()
file_path = os.path.dirname(output_path)
del_ins_type_size_chart=file_path+os.path.sep+"del_ins_type_size.jpg"
lenght_variant_file=file_path+os.path.sep+"lenght_variant.jpg"
del_ins_genotype_chart=file_path+os.path.sep+"del_ins_genotype.jpg"
genotype_chart1=file_path+os.path.sep+"inv_dup_genotype.jpg"
out_chart_2=file_path +os.path.sep+"dup_inv_type_size.jpg"
vcf_variables=vcf_number_variants(input_vcf_file)

INS_DUP=vcf_variables.INS+vcf_variables.DUP
sv_ranges["DEL"] = count_numbers_in_ranges(vcf_variables.DEL, ranges)
sv_ranges["INS"] = count_numbers_in_ranges(vcf_variables.INS, ranges)
sv_ranges["DUP"] = count_numbers_in_ranges(vcf_variables.DUP, ranges)
sv_ranges["INV"] = count_numbers_in_ranges(vcf_variables.INV, ranges)
sv_ranges["BND"] = count_numbers_in_ranges(vcf_variables.BND, ranges)
allele_frequency_chart_genrator(input_vcf_file, output_path)
sv_size_type_chart(sv_ranges["DEL"],sv_ranges["INS"],"DEL","INS",del_ins_type_size_chart)
sv_size_type_chart(sv_ranges["DUP"],sv_ranges["INV"],"DUP","INV",out_chart_2)
#genome_bar_chart((DEL_GENOTYPE,"DEL"),(INS_GENOTYPE,"INS"),del_ins_genotype_chart )
lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 10000, 50000000,1, 1, [10000,25000000,50000000],  0,"DEl>=10K\nbin_size=500K")
lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 1000, 10000, 2, 0, [1000,5000],  0,"1k<=DEl<=10K\nbin_size=90")
lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 50, 1000, 3, 0,[50,500],  0,"50<=DEl<=1K\nbin_size=10")
lenght_var_count_chart(lenght_variant_file,0,vcf_variables.INS, 50, 1000, 4, 0,[500,1000],  0,"50<=INS<=1K\nbin_size=10")
lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 1000, 10000, 5, 0,[5000,10000],  0,"1k<=INS/DUP<=10K\nbin_size=90")
lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 0, 50000000,6, 0, [25000000,50000000],  0,"INS/DUP>=10K\nbin_size=500K")
