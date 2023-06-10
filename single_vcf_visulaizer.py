# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC - Farhang Jaryani
"""
from arg_parser import argument_parser
from functions_variables_multi import *
from functions_variables_single import *
import os
sv_ranges = {
    "DEL": None,
    "INS": None,
    "DUP": None,
    "INV": None,
    "BND": None,
}
labels=["0/0","0/1","1/1"]

def single_visulaizer(input_vcf_file,output_path):
    #input_vcf_file,output_path=argument_parser()

    del_ins_type_size_chart=os.path.join(output_path,"del_ins_type_size.jpg")
    lenght_variant_file=os.path.join(output_path,"lenght_variant.jpg")
    genome_chart_del_ins=os.path.join(output_path,"del_ins_genotype.jpg")
    genotype_chart_inv_dup=os.path.join(output_path,"inv_dup_genotype.jpg")
    out_chart_2=os.path.join(output_path ,"dup_inv_type_size.jpg")
    vcf_variables=vcf_number_variants(input_vcf_file)

    INS_DUP=vcf_variables.INS+vcf_variables.DUP
    sv_ranges["DEL"] = count_numbers_in_ranges(vcf_variables.DEL, ranges)
    sv_ranges["INS"] = count_numbers_in_ranges(vcf_variables.INS, ranges)
    sv_ranges["DUP"] = count_numbers_in_ranges(vcf_variables.DUP, ranges)
    sv_ranges["INV"] = count_numbers_in_ranges(vcf_variables.INV, ranges)
    sv_ranges["BND"] = count_numbers_in_ranges(vcf_variables.BND, ranges)

    Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_path)
    Genome_chart_data_generator.allele_frequency_chart_generator()

    genome_bar_chart(genotype_chart_inv_dup,labels,GenomeChartData(vcf_variables.INV_GENOTYPE,"INV"),GenomeChartData(vcf_variables.DUP_GENOTYPE,"DUP"))

    Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_path)
    sv_size_type_chart(sv_ranges["DEL"],sv_ranges["INS"],"DEL","INS",del_ins_type_size_chart)
    sv_size_type_chart(sv_ranges["DUP"],sv_ranges["INV"],"DUP","INV",out_chart_2)
    #genome_bar_chart((DEL_GENOTYPE,"DEL"),(INS_GENOTYPE,"INS"),del_ins_genotype_chart )
    genome_bar_chart(genome_chart_del_ins,labels,GenomeChartData(vcf_variables.DEL_GENOTYPE,"DEl"),GenomeChartData(vcf_variables.INS_GENOTYPE,"INS"))


    lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 10000, 50000000,1, 1, [10000,25000000,50000000],  0,"DEl>=10K\nbin_size=500K")
    lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 1000, 10000, 2, 0, [1000,5000],  0,"1k<=DEl<=10K\nbin_size=90")
    lenght_var_count_chart(lenght_variant_file,1,vcf_variables.DEL, 50, 1000, 3, 0,[50,500],  0,"50<=DEl<=1K\nbin_size=10")
    lenght_var_count_chart(lenght_variant_file,0,vcf_variables.INS, 50, 1000, 4, 0,[500,1000],  0,"50<=INS<=1K\nbin_size=10")
    lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 1000, 10000, 5, 0,[5000,10000],  0,"1k<=INS/DUP<=10K\nbin_size=90")
    lenght_var_count_chart(lenght_variant_file,0,INS_DUP, 0, 50000000,6, 0, [25000000,50000000],  0,"INS/DUP>=10K\nbin_size=500K")

if __name__ == "__main__":
    main()
