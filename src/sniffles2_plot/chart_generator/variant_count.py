# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os

import matplotlib.pyplot as plt
import pandas as pd

from sniffles2_plot.helper.io_class import FileIO
from sniffles2_plot.parser.vcf_line_parser import VCFLineSVPopulation


class VariantCount(FileIO):
    """generate plots for variant count all SVs"""

    def variant_count_chart_generator(self):
        """generate the allele frequency plots"""
        self.del_count = 0
        self.dup_count = 0
        self.ins_count = 0
        self.inv_count = 0
        self.bnd_count = 0
        self.cnv_count = 0
        self.other_count = 0
        test_sv = []
        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                obj = VCFLineSVPopulation(line.strip())
                if obj.ERROR:
                    continue
                sv_type = obj.SVTYPE
                samples_AF = obj.samples_AF
                test_sv.append(sv_type)

                if sv_type == "DEL":
                    self.del_count += 1
                elif sv_type == "INS":
                    self.ins_count += 1
                elif sv_type == "INV":
                    self.inv_count += 1
                elif sv_type == "DUP":
                    self.dup_count += 1
                elif sv_type == "BND":
                    self.bnd_count += 1
                elif sv_type == "CNV":
                    self.cnv_count += 1
                else:
                    self.other_count += 1
        x_labels = ["DEL", "DUP", "INS", "INV", "BND", "CNV", "OTHER"]

        frequencies = [
            self.del_count,
            self.dup_count,
            self.ins_count,
            self.inv_count,
            self.bnd_count,
            self.cnv_count,
            self.other_count,
        ]
        # Create the bar plot

        colors = ["blue", "green", "orange", "red", "purple", "yellow", "cyan"]
        plt.bar(x_labels, frequencies, color=colors, label=x_labels)
        for i, value in enumerate(frequencies):
            plt.text(i, value, str(value), ha="center", va="bottom")
        # Add labels and title
        plt.xlabel("SV Types")
        plt.yscale("log")
        plt.ylabel("SV count")
        plt.title(f"Variant call (all SV)\n n={sum(frequencies)}")

        # Show the plot
        plt.savefig(self.output_file("variant_count.jpg"), dpi=800)
        plt.close()
