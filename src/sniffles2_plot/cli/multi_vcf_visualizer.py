# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC Farhang Jaryani
"""

from sniffles2_plot.chart_generator import (
    GenomeChartDataGenerator,
    SizeDistribution,
    Sv_sites_per_genome,
    VariantCount,
)


def generate_multi_vcf_charts(input_vcf_file, output_path):
    """Multi VCF plot generator"""
    SV_site = Sv_sites_per_genome(input_vcf_file, output_path)
    Genome_chart_data_generator = GenomeChartDataGenerator(input_vcf_file, output_path)
    S_D_obj = SizeDistribution(input_vcf_file, output_path)
    V_C_obj = VariantCount(input_vcf_file, output_path)

    try:
        SV_site.sv_sites_per_genome()
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        V_C_obj.variant_count_chart_generator()
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        Genome_chart_data_generator.allele_frequency_chart_generator()
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        Genome_chart_data_generator.samples_sv_numbers()
    except Exception as e:
        print(f"An error occurred: {e}")

    try:
        Genome_chart_data_generator.heat_map_generator()
    except Exception as e:
        print(f"An error occurred: {e}")

