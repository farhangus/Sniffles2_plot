# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os
import collections
import seaborn as sns
import pandas as pd
from upsetplot import plot
from upsetplot import from_memberships
import matplotlib.pyplot as plt
from vcf_line_parser import VCFLineSVPopulation
from DataClasses import *
from functions_variables_single import *

def sample_to_matrix(sample_names,samples):

    tmp_dict=count_frequency(samples)
    lenght=len(sample_names)
    string1=lenght*"0"
    tmp_matrix=[]
    for i in range(lenght):
        for j in range(lenght):
            tmp_string=list(string1)
            tmp_string[i]="1"
            tmp_string[j]="1"
            tmp_sum=0
            for key,value in tmp_dict.items():
                if key[i]=="1" and key[j]=="1":
                    tmp_sum+=value
            tmp_matrix.append(tmp_sum)
    two_dim_matrix = np.reshape(tmp_matrix, (lenght,lenght))
    return two_dim_matrix

def count_frequency(elements):
    frequency_dict = {}
    for element in elements:
        if element in frequency_dict:
            frequency_dict[element] += 1
        else:
            frequency_dict[element] = 1
    return frequency_dict

def samples_SV_counter(input_file_name, output_results):
    """counting the number of SVs"""
    with open(input_file_name, 'r',encoding="utf-8") as f:
        c = collections.Counter(text.strip() for text in f)
    with open(output_results, 'w',encoding="utf-8") as f:
        f.write("\n".join(
            f"{val} {key}" for key, val in c.most_common()))

class GenomeChartDataGenerator:
    """generate plots for multi VCF files"""
    def __init__(self, input_file_path, output_directory):
        self.input_file_path = input_file_path
        self.output_directory = output_directory

    def output_file(self, filename):
        """returnthe output file name and filepath"""
        return os.path.join(self.output_directory,filename)

    def allele_frequency_chart_generator(self):
        """generate the allele frequency plots"""
        sample_names = []
        with open(self.input_file_path, "r",encoding="utf-8") as f:
            for line in f:
                if line.startswith("#C"):
                    sample_names = [l.strip() for l in line.split('\t')[9:]]
                    break
        if not sample_names:
            sample_names = ['sample_1']

        DEL_LIST = []
        INS_LIST = []
        INV_LIST = []
        DUP_LIST = []

        with open(self.input_file_path, "r",encoding="utf-8") as f:
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

        NUMBER_SAMPLES = len(DEL_LIST[0][:1000])

        bin_size = 0.04
        num_bins = int(1 / bin_size)

        for i in range(NUMBER_SAMPLES):
            tmp_list_name_del = separate_lists(DEL_LIST, i)
            tmp_list_name_ins = separate_lists(INS_LIST, i)
            tmp_list_name_inv = separate_lists(INV_LIST, i)
            tmp_list_name_dup = separate_lists(DUP_LIST, i)

            plt.hist([tmp_list_name_del, tmp_list_name_ins, tmp_list_name_inv, tmp_list_name_dup],
                     bins=num_bins, range=(0, 1), label=['DEL', 'INS', 'INV', 'DUP'],
                     alpha=0.7, edgecolor='black')

            plt.yscale('log')
            plt.xlabel(sample_names[i])
            plt.ylabel('Count (log)')
            plt.title('Variant Frequency Spectrum')
            plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1), prop={'size': 5})

            plt.savefig(self.output_file(sample_names[i]), dpi=800)
            plt.close()

    def samples_sv_numbers(self):
        """generate the upset plots for SVs per each sample"""
        sample_names=[]
        with open(self.output_file("tmp.txt"), "w",encoding="utf-8") as f_sv_samples:
            with open(self.input_file_path, "r",encoding="utf-8") as f:
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
                            f_sv_samples.write(obj.SUPP_VEC+"\n")
        if not sample_names:
            sample_names = ['sample_1']

        samples_SV_counter(
            self.output_file("tmp.txt"),
            self.output_file("sv_sample_results.txt"))
        data = []
        data_set = []
        with open(self.output_file("sv_sample_results.txt"), "r",encoding="utf-8") as f:
            for line in f:
                elements = line.split()
                tmp_line = elements[-1].strip("\n")
                data_set.append(tmp_line)
                data.append(int(elements[0]))

        data_list = []
        for item in data_set:
            if item != "0":
                tmp_list = [
                    sample_names[i]
                    for i, item_i in enumerate(item)
                    if item_i != "0"
                ]
                data_list.append(tmp_list)
            else:
                data_list.append([])

        example = from_memberships(data_list, data=data)
        upset = plot(example, facecolor="black", other_dots_color=.4, shading_color=.1)

        for patch in upset['intersections'].patches:
            upset['intersections'].annotate(
                text=patch.get_height(),
                xy=(patch.get_x() + patch.get_width() / 2, patch.get_height()),
                ha='center',
                va='bottom',
                rotation='vertical',
                fontsize=6,
                xytext=(0, 5),
                textcoords='offset points')

        upset['intersections'].bar_color = 'blue'
        upset['intersections'].bar_alpha = 0.7
#        os.remove(self.output_file("tmp.txt"))
#        os.remove(self.output_file("sv_sample_results.txt"))
        plt.savefig(self.output_file("sample_upset.png"), dpi=800, edgecolor="white")

    def heat_map_generator(self):
        print("heat ap generator")
        sample_names=[]
        samples_del=[]
        samples_ins=[]
        with open(self.output_file("tmp.txt"), "w",encoding="utf-8") as f_sv_samples:
            with open(self.input_file_path, "r",encoding="utf-8") as f:
                for line in f:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#C"):
                        sample_names = [l.strip() for l in line.split('\t')[9:]]
                    else:
                        obj = VCFLineSVPopulation(line)
                        if obj.ERROR:
                            continue
                        if obj.FILTER=="PASS" and obj.SVTYPE=="DEL":
                            f_sv_samples.write(f"{obj.SUPP_VEC}\n")
                            samples_del.append(obj.SUPP_VEC)
                        elif obj.FILTER=="PASS" and obj.SVTYPE=="INS":
                            samples_ins.append(obj.SUPP_VEC)

        del_matrix=sample_to_matrix(sample_names,samples_del)
        ins_matrix=sample_to_matrix(sample_names,samples_ins)

        print(del_matrix)
        print(ins_matrix)
        combined_matrix=np.tril(del_matrix) + np.tril(ins_matrix, -1).transpose()
        print(combined_matrix)
        data = np.array(combined_matrix)
        # # mask = np.triu(np.ones_like(data))
        # df = pd.DataFrame(combined_matrix, columns=sample_names)
        df_cm = pd.DataFrame(combined_matrix, index = sample_names,
                  columns = sample_names)
        plt.figure(figsize=(10,10))

        sns.heatmap(df_cm, cmap='PuOr')


        # # Create a heatmap
        # sns.heatmap(data, annot=True, fmt="d", cmap="YlGnBu")

        # # Set x and y axis labels
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        # plt.yticks(np.arange(9), range(9))

        # # Display the heatmap
        plt.savefig("hetamap.jpg",dpi=800)
