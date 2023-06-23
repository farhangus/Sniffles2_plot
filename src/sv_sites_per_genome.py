# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from vcf_line_parser import VCFLineSVPopulation


class Sv_sites_per_genome:
    """generate plots for variant count all SVs"""

    def __init__(self, input_file_path, output_directory):
        self.input_file_path = input_file_path
        self.output_directory = output_directory

    def output_file(self, filename):
        """returnthe output file name and filepath"""
        return os.path.join(self.output_directory, filename)

    def sv_sites_per_genome(self):
        """generate the allele frequency plots"""
        self.del_dict = 0
        self.dup_dict = 0
        self.ins_dict = 0
        self.inv_dict = 0
        self.bnd_dict = 0
        self.cnv_dict = 0
        self.other_dict = 0
        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#C"):
                    sample_names = [l.strip() for l in line.split('\t')[9:]]
                else:
                    obj = VCFLineSVPopulation(line)
                    if obj.ERROR:
                        continue
                    if obj.FILTER=="PASS":
                        print(obj.SVTYPE,obj.samples_AF)
                        # if obj.sv_type == "DEL":
                        #     self.del_dict []
                        # elif sv_type == "INS":
                        #     self.ins_dict += 1
                        # elif sv_type == "INV":
                        #     self.inv_dict += 1
                        # elif sv_type == "DUP":
                        #     self.dup_dict += 1
                        # elif sv_type == "BND":
                        #     self.bnd_dict += 1
                        # elif sv_type == "CNV":
                        #     self.cnv_dict += 1
                        # else:
                        #     self.other_dict += 1
        # x_labels = ["DEL", "DUP", "INS", "INV", "BND", "CNV", "OTHER"]

        # frequencies = [
        #     self.del_count,
        #     self.dup_count,
        #     self.ins_count,
        #     self.inv_count,
        #     self.bnd_count,
        #     self.cnv_count,
        #     self.other_count,
        # ]
        # # Create the bar plot

        # print(x_labels)
        # colors = ["blue", "green", "orange", "red", "purple", "yellow", "cyan"]
        # print(type(x_labels))
        # plt.bar(x_labels, frequencies, color=colors,label = x_labels)
        # for i, value in enumerate(frequencies):
        #     plt.text(i, value, str(value), ha="center", va="bottom")
        # # Add labels and title
        # plt.xlabel("SV Types")
        # plt.yscale("log")
        # plt.ylabel("SV count")
        # plt.title(f"Variant call (all SV)\n n={sum(frequencies)}")

        # # Show the plot
        # plt.savefig(self.output_file("variant_count.jpg"),dpi=800)
        # plt.close()