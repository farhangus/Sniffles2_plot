# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC
"""
from arg_parser import argument_parser
from functions_variables import *
import os
input_vcf_file,output_chrt=argument_parser()
file_path = os.path.dirname(output_chrt)
lenght_variant_file=file_path+os.path.sep+"lenght_variant.jpg"
genotype_chart=file_path+os.path.sep+"genotype.jpg"
genotype_chart1=file_path+os.path.sep+"genotype_inv_dup.jpg"
out_chart_2=file_path +os.path.sep+"sv_size_type_dup_inv.jpg"
DEL,DEL_GENOTYPE,INS,INS_GENOTYPE,DUP,DUP_GENOTYPE,INV,INV_GENOTYPE,BND,BND_GENOTYPE=vcf_number_variants(input_vcf_file)
genome_bar_chart((DUP_GENOTYPE,"DUP"),(INV_GENOTYPE,"INV"),genotype_chart1 )

INS_DUP=INS+DUP
Result_DEL = count_numbers_in_ranges(DEL, ranges)
Result_INS = count_numbers_in_ranges(INS, ranges)
Result_DUP = count_numbers_in_ranges(DUP, ranges)
Result_INV = count_numbers_in_ranges(INV, ranges)
Result_BND = count_numbers_in_ranges(BND, ranges)
sv_size_type_chart(Result_DEL,Result_INS,"DEL","INS",output_chrt)
sv_size_type_chart(Result_DUP,Result_INV,"DUP","INV",out_chart_2)
genome_bar_chart((DEL_GENOTYPE,"DEL"),(INS_GENOTYPE,"INS"),genotype_chart )

lenght_var_count_chart(lenght_variant_file,1,DEL, 50, 1000, 3, 0,[50,500],  0,"50<=DEl<=1K")
lenght_var_count_chart(lenght_variant_file,1,DEL, 1000, 10000, 2, 0, [1000,5000],  0,"1k<=DEl<=10K")
lenght_var_count_chart(lenght_variant_file,1,DEL, 10000, 50000000,1, 1, [10000,25000000,50000000],  0,"DEl>=10K")
lenght_var_count_chart(lenght_variant_file,0,INS, 50, 1000, 4, 0,[500,1000],  0,"50<=INS<=1K")
lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 1000, 10000, 5, 0,[5000,10000],  0,"1k<=INS/DUP<=10K")
lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 0, 50000000,6, 0, [25000000,50000000],  0,"INS/DUP>=10K")
