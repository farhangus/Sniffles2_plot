# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC - Farhang Jaryani
"""
import os

from sniffles2_plot.chart_generator import (
    GenomeChartData,
    GenomeChartDataGenerator,
    VariantCount,
    count_numbers_in_ranges,
    genome_bar_chart,
    length_var_count_chart,
    ranges,
    sv_size_type_chart,
    vcf_number_variants,
)

sv_ranges = {
    "DEL": None,
    "INS": None,
    "DUP": None,
    "INV": None,
    "BND": None,
}
labels = ["0/0", "0/1", "1/1"]
phased_labels=["0|0","0|1","1|1"]

def single_visulaizer(input_vcf_file, output_path):
    """generate the plots for single vcf files"""
    V_C_obj = VariantCount(input_vcf_file, output_path)
    V_C_obj.variant_count_chart_generator()

    del_ins_type_size_chart = os.path.join(output_path, "del_ins_type_size.jpg")
    length_variant_file = os.path.join(output_path, "length_variant.jpg")
    genome_chart_del_ins = os.path.join(output_path, "del_ins_genotype.jpg")
    genotype_chart_inv_dup = os.path.join(output_path, "inv_dup_genotype.jpg")
    genome_chart_phased_del_ins = os.path.join(output_path, "phased_del_ins_genotype.jpg")
    genotype_chart_phased_inv_dup = os.path.join(output_path, "phased_inv_dup_genotype.jpg")
    out_chart_2 = os.path.join(output_path, "dup_inv_type_size.jpg")
    vcf_variables = vcf_number_variants(input_vcf_file)

    INS_DUP = vcf_variables.INS + vcf_variables.DUP
    sv_ranges["DEL"] = count_numbers_in_ranges(vcf_variables.DEL, ranges)
    sv_ranges["INS"] = count_numbers_in_ranges(vcf_variables.INS, ranges)
    sv_ranges["DUP"] = count_numbers_in_ranges(vcf_variables.DUP, ranges)
    sv_ranges["INV"] = count_numbers_in_ranges(vcf_variables.INV, ranges)
    sv_ranges["BND"] = count_numbers_in_ranges(vcf_variables.BND, ranges)

    Genome_chart_data_generator = GenomeChartDataGenerator(input_vcf_file, output_path)
    Genome_chart_data_generator.allele_frequency_chart_generator()

    genome_bar_chart(
        genotype_chart_inv_dup,
        labels,
        GenomeChartData(vcf_variables.INV_GENOTYPE, "INV"),
        GenomeChartData(vcf_variables.DUP_GENOTYPE, "DUP"),
    )


    Genome_chart_data_generator = GenomeChartDataGenerator(input_vcf_file, output_path)
    sv_size_type_chart(
        sv_ranges["DEL"], sv_ranges["INS"], "DEL", "INS", del_ins_type_size_chart
    )
    sv_size_type_chart(sv_ranges["DUP"], sv_ranges["INV"], "DUP", "INV", out_chart_2)
    # genome_bar_chart((DEL_GENOTYPE,"DEL"),(INS_GENOTYPE,"INS"),del_ins_genotype_chart )
    genome_bar_chart(
        genome_chart_del_ins,
        labels,
        GenomeChartData(vcf_variables.DEL_GENOTYPE, "DEl"),
        GenomeChartData(vcf_variables.INS_GENOTYPE, "INS"),
    )
    if vcf_variables.HAS_PHASED:
        genome_bar_chart(
            genome_chart_phased_del_ins,
            phased_labels,
            GenomeChartData(vcf_variables.PHASED_DEL_GENOTYPE, "PHASED_DEl"),
            GenomeChartData(vcf_variables.PHASED_INS_GENOTYPE, "PHASED_INS"),
        )
        genome_bar_chart(
            genotype_chart_phased_inv_dup,
            phased_labels,
            GenomeChartData(vcf_variables.PHASED_INV_GENOTYPE, "PHASED_INV"),
            GenomeChartData(vcf_variables.PHASED_DUP_GENOTYPE, "PHASED_DUP"),
        )
        
    length_var_count_chart(
        length_variant_file,
        1,
        vcf_variables.DEL,
        10000,
        50000000,
        1,
        1,
        [25000000, 50000000],
        0,
        "DEl>=10K\nbin_size=500K",
    )
    length_var_count_chart(
        length_variant_file,
        1,
        vcf_variables.DEL,
        1000,
        10000,
        2,
        0,
        [5000, 10000],
        0,
        "1k<=DEl<=10K\nbin_size=90",
    )
    length_var_count_chart(
        length_variant_file,
        1,
        vcf_variables.DEL,
        50,
        1000,
        3,
        0,
        [0,500, 1000],
        0,
        "50<=DEl<=1K\nbin_size=10",
    )
    length_var_count_chart(
        length_variant_file,
        0,
        vcf_variables.INS,
        50,
        1000,
        4,
        0,
        [500, 1000],
        0,
        "50<=INS<=1K\nbin_size=10",
    )
    length_var_count_chart(
        length_variant_file,
        0,
        INS_DUP,
        1000,
        10000,
        5,
        0,
        [5000, 10000],
        0,
        "1k<=INS/DUP<=10K\nbin_size=90",
    )
    length_var_count_chart(
        length_variant_file,
        0,
        INS_DUP,
        10000,
        50000000,
        6,
        0,
        [25000000, 50000000],
        0,
        "INS/DUP>=10K\nbin_size=500K",
    )
