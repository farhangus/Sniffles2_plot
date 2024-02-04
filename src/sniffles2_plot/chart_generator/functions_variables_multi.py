# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import collections
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from upsetplot import from_memberships, plot

from sniffles2_plot.chart_generator.functions_variables_single import *
from sniffles2_plot.parser.vcf_line_parser import VCFLineSVPopulation
from sniffles2_plot.schemas import *

sys.path.append("src/")
from sniffles2_plot.helper.io_class import FileIO


def sample_to_matrix(sample_names, samples):
    tmp_dict = count_frequency(samples)
    length = len(sample_names)
    string1 = length * "0"
    tmp_matrix = []
    for i in range(length):
        for j in range(length):
            tmp_string = list(string1)
            tmp_string[i] = "1"
            tmp_string[j] = "1"
            tmp_sum = 0
            for key, value in tmp_dict.items():
                if key[i] == "1" and key[j] == "1":
                    tmp_sum += value
            tmp_matrix.append(tmp_sum)
    two_dim_matrix = np.reshape(tmp_matrix, (length, length))
    return two_dim_matrix


def count_frequency(elements):
    frequency_dict = {}
    for element in elements:
        if element in frequency_dict:
            frequency_dict[element] += 1
        else:
            frequency_dict[element] = 1
    return frequency_dict


def samples_SV_counter(output_results, samples_intersections=()):
    """counting the number of SVs"""
    c_ = collections.Counter(samples_intersections)
    with open(output_results, "w", encoding="utf-8") as f:
        f.write("\n".join(f"{val} {key}" for key, val in c_.most_common()))


class GenomeChartDataGenerator(FileIO):
    """generate plots for multi VCF files"""

    def allele_frequency_chart_generator(self):
        """generate the allele frequency plots"""
        sample_names = []
        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#C"):
                 #   sample_names = [l.strip() + str(index+1) for index, l in enumerate(line.split("\t")[9:])]
                    sample_names = line.split("\t")[9:]
                    break
        if not sample_names:
            sample_names = ["sample_1"]
        DEL_LIST = []
        INS_LIST = []
        INV_LIST = []
        DUP_LIST = []

        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                obj = VCFLineSVPopulation(line.strip())
                if obj.ERROR:
                    continue
                sv_type = obj.SVTYPE
                samples_AF = obj.samples_AF

                if sv_type == "DEL":
                    DEL_LIST.append(samples_AF)
                elif sv_type == "INS":
                    INS_LIST.append(samples_AF)
                elif sv_type == "INV":
                    INV_LIST.append(samples_AF)
                elif sv_type == "DUP":
                    DUP_LIST.append(samples_AF)
        NUMBER_SAMPLES = len(sample_names)
        bin_size = 0.04
        num_bins = int(1 / bin_size)
        for i in range(NUMBER_SAMPLES):
            tmp_list_name_del = separate_lists(DEL_LIST, i)
            tmp_list_name_ins = separate_lists(INS_LIST, i)
            tmp_list_name_inv = separate_lists(INV_LIST, i)
            tmp_list_name_dup = separate_lists(DUP_LIST, i)
            plt.hist(
                [
                    tmp_list_name_del,
                    tmp_list_name_ins,
                    tmp_list_name_inv,
                    tmp_list_name_dup,
                ],
                bins=num_bins,
                range=(0, 1),
                label=["DEL", "INS", "INV", "DUP"],
                alpha=0.7,
                edgecolor="black",
            )

            plt.yscale("log")
            plt.xlabel(sample_names[i])
            plt.ylabel("Count (log)")
            plt.title("Variant Frequency Spectrum")
            plt.legend(loc="upper right", bbox_to_anchor=(1.1, 1.1), prop={"size": 5})
            plt.savefig(self.output_file(sample_names[i].replace(".", "_")), dpi=800)
            plt.close()

    def samples_sv_numbers(self):
        """generate the upset plots for SVs per each sample"""
        sample_names = []
        samples_intersections=[]
        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#C"):
                    sample_names = [l.strip() for index, l in enumerate(line.split("\t")[9:])]
                else:
                    obj = VCFLineSVPopulation(line)
                    if obj.ERROR:
                        continue
                    if obj.FILTER == "PASS":
                        samples_intersections.append(obj.SUPP_VEC)
        if not sample_names:
            sample_names = ["sample_1"]
        samples_SV_counter(
            self.output_file("sv_sample_results.txt"),samples_intersections
        )
        data = []
        data_set = []

        with open(
            self.output_file("sv_sample_results.txt"), "r", encoding="utf-8"
        ) as f:
            for line in f:
                elements = line.split()
                tmp_line = elements[-1].strip("\n")
                data_set.append(tmp_line)
                data.append(int(elements[0]))
        data_list = []
        for item in data_set:
            if item != "0":
                tmp_list = [
                    sample_names[i] for i, item_i in enumerate(item) if item_i != "0"
                ]

               # tmp_list = sample_names
                data_list.append(tmp_list)
            else:
                data_list.append([])
        
        if len(data_list) < 2:
            return 
        
        example = from_memberships(data_list, data=data)
        
        upset = plot(
            example, facecolor="black", other_dots_color=0.4, shading_color=0.1
        )

        for patch in upset["intersections"].patches:
            upset["intersections"].annotate(
                text=patch.get_height(),
                xy=(patch.get_x() + patch.get_width() / 2, patch.get_height()),
                ha="center",
                va="bottom",
                rotation="vertical",
                fontsize=6,
                xytext=(0, 5),
                textcoords="offset points",
            )
        upset["intersections"].bar_color = "blue"
        upset["intersections"].bar_alpha = 0.7
        os.remove(self.output_file("sv_sample_results.txt"))
        plt.savefig(self.output_file("sample_upset.jpg"), dpi=400, edgecolor="white")
        plt.close()

    def heat_map_generator(self):
        sample_names = []
        samples_del = []
        samples_ins = []
        with open(self.input_file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#C"):
                  #  sample_names = [l.strip() + str (index+1) for index, l in enumerate(line.split("\t")[9:])]
                    sample_names = [l.strip()  for index, l in enumerate(line.split("\t")[9:])]
                else:
                    obj = VCFLineSVPopulation(line)
                    if obj.ERROR:
                        continue
                    if obj.FILTER == "PASS" and obj.SVTYPE == "DEL":
                        samples_del.append(obj.SUPP_VEC)
                    elif obj.FILTER == "PASS" and obj.SVTYPE == "INS":
                        samples_ins.append(obj.SUPP_VEC)
        # counter = sum(1 for i in samples_ins if i[0] == "1" and i[-1] == "1")
        # print(counter)
        del_matrix = sample_to_matrix(sample_names, samples_del)
        ins_matrix = sample_to_matrix(sample_names, samples_ins)
        combined_matrix = np.tril(del_matrix) + np.tril(ins_matrix, -1).transpose()

        data = np.array(combined_matrix)
        np.fill_diagonal(combined_matrix, 0)

        # # mask = np.triu(np.ones_like(data))
        # df = pd.DataFrame(combined_matrix, columns=sample_names)
        df_cm = pd.DataFrame(combined_matrix, index=sample_names, columns=sample_names)
        plt.figure(figsize=(15,15))
        sns.heatmap(df_cm, cmap="PuOr", annot=False, fmt=".0f")
        plt.yticks(rotation=0)
        plt.xticks(rotation=90)
        # # Create a heatmap
        # sns.heatmap(data, annot=True, fmt="d", cmap="YlGnBu")
        # # Set x and y axis labels
        plt.title("right Insertion & left Deletion")
        plt.xlabel("sample_names")
        plt.ylabel("sample_names")
        n = len(sample_names)
        line_color = (0, 0, 1, 0.1)  # RGBA color tuple (R=0, G=0, B=1, alpha=0.5)
        plt.plot([0, n], [0, n], color=line_color, linewidth=2)
        # plt.yticks(np.arange(9), range(9))
        plt.tight_layout()
        # # Display the heatmap
        plt.savefig(self.output_file("heatmap.jpg"), dpi=400)
        plt.close()
