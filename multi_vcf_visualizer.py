# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""
import sys
sys.path.append('src/')  
from  variant_count import *
from functions_variables_multi import GenomeChartDataGenerator
from size_distribution import *

def multi_visulaizer(input_vcf_file,output_path):
    """Multi VCF plot generator"""
    V_C_obj=VariantCount(input_vcf_file, output_path)
    V_C_obj.variant_count_chart_generator()
    S_D_obj=SizeDistribution(input_vcf_file, output_path)
    S_D_obj.generate_size_distribution_plot()
    # Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_chrt)
    # Genome_chart_data_generator.allele_frequency_chart_generator()
    # Genome_chart_data_generator.samples_sv_numbers()
    # Genome_chart_data_generator.heat_map_generator()
    
