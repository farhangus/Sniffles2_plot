# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from sniffles2_plot.parser.vcf_line_parser import VCFLineSV
from sniffles2_plot.schemas import *

SAMPLES = 0
from sniffles2_plot.helper.io_class import FileIO


def list_percentage(list_):
    total_sum = sum(list_)
    modified_list = [(x / total_sum) * 100 if x != 0 else 0 for x in list_]
    return modified_list


@dataclass
class VcfVariables:
    del_size: list
    dup_size: list
    cnv_size: list
    ins_size: list
    inv_size: list
    cpx_size: list
    all_size: list

    @classmethod
    def new(cls):
        return cls([], [], [], [], [], [], [])


class SizeDistribution(FileIO):
    """generate size distribution plots"""

    def __init__(self, *args):
        super().__init__(*args)
        self.__private_counter = 0

    def variants_couns(self):
        """return the number of each variant seprately"""
        global SAMPLES
        with open(self.input_file_path, "r") as f:
            lines = f.readlines()
            vcf_variables = VcfVariables.new()
            for line in lines:
                if line[0] != "#":
                    obj = VCFLineSV(line)
                    if obj.ERROR:
                        continue
                    SAMPLES += 1
                    vcf_variables.all_size.append(abs(obj.SVLEN))
                    if obj.SVTYPE == "DEL":
                        vcf_variables.del_size.append(abs(obj.SVLEN))
                    elif obj.SVTYPE == "INS":
                        vcf_variables.ins_size.append(obj.SVLEN)
                    elif obj.SVTYPE == "INV":
                        vcf_variables.inv_size.append(obj.SVLEN)
                    elif obj.SVTYPE == "DUP":
                        vcf_variables.dup_size.append(obj.SVLEN)
                    elif obj.SVTYPE == "CNV":
                        vcf_variables.cnv_size.append(obj.SVLEN)
                    elif obj.SVTYPE == "CPX":
                        vcf_variables.cpx_size.append(obj.SVLEN)
        return vcf_variables

    def generate_size_distribution_plot(self):
        tmp = self.variants_couns()
        all_size = result = [int(x / 6) for x in tmp.all_size]
        counts, bins = np.histogram(all_size, bins=np.logspace(0, 6, 50))
        bin_centers = (bins[:-1] + bins[1:]) / 2
        counts = list_percentage(counts)
        plt.plot(bin_centers, counts, label="ALL", color="black", linewidth=1.5)
        plt.scatter(bin_centers, counts, color="black", s=5)

        counts, bins = np.histogram(tmp.del_size, bins=np.logspace(0, 6, 50))
        counts = list_percentage(counts)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        plt.plot(bin_centers, counts, label="DEL", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        counts, bins = np.histogram(tmp.ins_size, bins=np.logspace(0, 6, 50))
        counts = list_percentage(counts)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        plt.plot(bin_centers, counts, label="INS", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        counts, bins = np.histogram(tmp.cpx_size, bins=np.logspace(0, 6, 50))
        bin_centers = (bins[:-1] + bins[1:]) / 2
        counts = list_percentage(counts)
        plt.plot(bin_centers, counts, label="CPX", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        counts, bins = np.histogram(tmp.inv_size, bins=np.logspace(0, 6, 50))
        counts = list_percentage(counts)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        plt.plot(bin_centers, counts, label="INV", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        counts, bins = np.histogram(tmp.dup_size, bins=np.logspace(0, 6, 50))
        counts = list_percentage(counts)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        plt.plot(bin_centers, counts, label="DUP", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        counts, bins = np.histogram(tmp.cnv_size, bins=np.logspace(0, 6, 50))
        counts = list_percentage(counts)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        plt.plot(bin_centers, counts, label="CNV", linewidth=0.5)
        plt.scatter(bin_centers, counts, s=5)

        # Set x-axis scale to logarithmic
        plt.xscale("log")
        plt.legend()
        text = f"N={SAMPLES}"
        # Set labels and title
        plt.xlabel("SIZE")

        plt.ylabel("Fraction of SV")
        plt.title(f"Size Distribution (ALL SV) \n {text} ")

        # Display the plot
        plt.savefig(self.output_file("size_distribution.jpg"), dpi=800)
        plt.close()
