# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import matplotlib.pyplot as plt
import numpy as np
from math import ceil

from sniffles2_plot.parser.vcf_line_parser import VCFLineSV
from sniffles2_plot.schemas import *

ranges = [
    (50, 100),
    (100, 200),
    (200, 300),
    (300, 400),
    (400, 600),
    (600, 800),
    (800, 1000),
    (1000, 2500),
    (2500, 5000),
    (5000, 10000000),
]
x_labels = [
    "50-100",
    "100-200",
    "200-300",
    "300-400",
    "400-600",
    "600-800",
    "800-1k",
    "1k-2.5k",
    "2.5k-5k",
    ">5k",
]


def separate_lists(lst, num):
    elements = [item[num] for item in lst]
    return elements


def count_numbers_in_ranges(numbers, ranges):
    counts = [0] * len(ranges)
    for number in numbers:
        for i, range_ in enumerate(ranges):
            if range_[0] <= number <= range_[1]:
                counts[i] += 1
                break
    return counts


def genome_bar_chart(output_file_path, labels, *bars: GenomeChartData):
    """genaret geome bar charts for (ins,del) and (dup,inv)"""
    x = np.arange(len(labels))
    # Define the width of each bar
    width = 0.2
    # Create the grouped bar plot
    plt.bar(x - width, bars[0].count(labels), width=width, label=bars[0].legend)
    plt.bar(x, bars[1].count(labels), width=width, label=bars[1].legend)
    # Add labels and title
    plt.xlabel("Genotype", size=6)
    plt.ylabel("Count", size=6)
    plt.title("Genotype_Frequency")
    # Set the x-axis tick labels
    plt.xticks(x, labels)
    # Add a legend
    plt.legend()
    plt.savefig(output_file_path, dpi=800)


def sv_size_type_chart(l1, l2, label_1, label_2, path_2_chart):
    """generate the SV size charts for each two variants(del,ins) and (inv,dup)"""
    fig, ax = plt.subplots()
    x = np.arange(len(l1))
    bar_width = 0.40
    plt.bar(x, l1, width=bar_width, label=label_1)
    plt.bar(x + bar_width, l2, width=bar_width, label=label_2)
    plt.xlabel("SizeBin")
    plt.ylabel("Count")
    plt.title("SV Size/Type Distribution")
    plt.xticks(x + bar_width / 5, x)
    ax.set_xticklabels(x_labels, rotation=45)
    plt.tight_layout()  # Adjust the padding to accommodate all labels
    plt.legend(title="svtype")
    plt.savefig(path_2_chart, dpi=800)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def length_var_count_chart(
    output_name,
    del_flag,
    v_s_list,
    min_value,
    max_value,
    chart_pos,
    ytick_flag,
    xtick_locs_list,
    plt_titl,
    subplt_title,
):
    """generate the comparision chart including 6 different charts for length variant"""
    plt.subplot(1, 6, chart_pos)
    plt.subplots_adjust(wspace=0.05, hspace=0.05, bottom=0.3)
    minimum = min_value
    maximum = max_value
    selected_values = [x for x in v_s_list if minimum <= x <= maximum]
    if len(selected_values) < 1:
        return
    bin_size =ceil((maximum - minimum)/100)
    # edge case1, len(selected values) < 2
    # len (set(selected_vales)) == 1
    # Create bins based on the bin size
    if bin_size == 0:
        bins = []
    else:
        bins = np.arange(minimum,maximum, bin_size)
    # Count the number of values that fall into each bin
    bin_counts, _ = np.histogram(selected_values, bins)
    # Display the bin plot
    plt.bar(bins[:-1], bin_counts, width=bin_size, align="edge")
    xtick_labels = [str(x) for x in xtick_locs_list]
    plt.xticks(xtick_locs_list, xtick_labels, rotation=90)
    if chart_pos == 1:
        plt.ylabel("Counts")
    if chart_pos == 3:
        plt.xlabel("\n\n        Lenght of Variants")
    plt.yscale("log")
    if ytick_flag:
        plt.yticks([1, 10, 100, 1000, 10000])
    else:
        plt.yticks([1, 10, 100, 1000, 10000])
        plt.yticks([])
        
    plt.title(subplt_title, size=5)
    if del_flag:
        plt.gca().invert_xaxis()
    plt.savefig(output_name, dpi=400)


def vcf_number_variants(input_vcf_file):
    """return the number of each variant seprately"""
    with open(input_vcf_file, "r") as f:
        lines = f.readlines()
        vcf_variables = VcfVariables.new()
        for line in lines:
            if line[0] != "#":
                obj = VCFLineSV(line)
                if obj.ERROR:
                    continue
                if obj.phased:
                    vcf_variables.HAS_PHASED=True
                if obj.SVTYPE == "DEL":
                    if obj.phased:
                        vcf_variables.PHASED_DEL.append(abs(obj.SVLEN))
                        vcf_variables.PHASED_DEL_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.PHASED_DEL_AF.append(obj.AF)
                    else:
                        vcf_variables.DEL.append(abs(obj.SVLEN))
                        vcf_variables.DEL_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.DEL_AF.append(obj.AF)
                if obj.SVTYPE == "INS":
                    if obj.phased:
                        vcf_variables.PHASED_INS.append(abs(obj.SVLEN))
                        vcf_variables.PHASED_INS_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.PHASED_INS_AF.append(obj.AF)
                    else:
                        vcf_variables.INS.append(abs(obj.SVLEN))
                        vcf_variables.INS_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.INS_AF.append(obj.AF)
                if obj.SVTYPE == "INV":
                    if obj.phased:
                        vcf_variables.PHASED_INV.append(abs(obj.SVLEN))
                        vcf_variables.PHASED_INV_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.PHASED_INV_AF.append(obj.AF)
                    else:
                        vcf_variables.INV.append(abs(obj.SVLEN))
                        vcf_variables.INV_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.INV_AF.append(obj.AF)
                if obj.SVTYPE == "DUP":
                    if obj.phased:
                        vcf_variables.PHASED_DUP.append(abs(obj.SVLEN))
                        vcf_variables.PHASED_DUP_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.PHASED_DUP_AF.append(obj.AF)
                    else:
                        vcf_variables.DUP.append(abs(obj.SVLEN))
                        vcf_variables.DUP_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.DUP_AF.append(obj.AF)
                if obj.SVTYPE == "BND":
                    if obj.phased:
                        vcf_variables.PHASED_BND.append(abs(obj.SVLEN))
                        vcf_variables.PHASED_BND_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.PHASED_BND_AF.append(obj.AF)
                    else:
                        vcf_variables.BND.append(abs(obj.SVLEN))
                        vcf_variables.BND_GENOTYPE.append(obj.GENOTYPE)
                        vcf_variables.BND_AF.append(obj.AF)

    return vcf_variables
