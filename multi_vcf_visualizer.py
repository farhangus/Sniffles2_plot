# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""
from functions_variables_multi import GenomeChartDataGenerator
def multi_visulaizer(input_vcf_file,output_chrt):
    """Multi VCF plot generator"""

    Genome_chart_data_generator=GenomeChartDataGenerator(input_vcf_file,output_chrt)
    Genome_chart_data_generator.allele_frequency_chart_generator()
    Genome_chart_data_generator.samples_sv_numbers()
    Genome_chart_data_generator.heat_map_generator()
    
