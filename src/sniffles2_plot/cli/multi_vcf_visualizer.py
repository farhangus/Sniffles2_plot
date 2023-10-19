# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""
import sys

sys.path.append("src/")
from sniffles2_plot.chart_generator import Sv_sites_per_genome,VariantCount,SizeDistribution,GenomeChartDataGenerator

def multi_visulaizer(input_vcf_file, output_path):
    """Multi VCF plot generator"""
    SV_site = Sv_sites_per_genome(input_vcf_file, output_path)
    SV_site.sv_sites_per_genome()
    V_C_obj = VariantCount(input_vcf_file, output_path)
    V_C_obj.variant_count_chart_generator()
    S_D_obj = SizeDistribution(input_vcf_file, output_path)
    #S_D_obj.generate_size_distribution_plot()
    Genome_chart_data_generator = GenomeChartDataGenerator(input_vcf_file, output_path)
    Genome_chart_data_generator.allele_frequency_chart_generator()
    Genome_chart_data_generator.samples_sv_numbers()
    Genome_chart_data_generator.heat_map_generator()
