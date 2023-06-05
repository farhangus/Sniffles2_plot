# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from vcf_line_parser import VCFLineSV
from vcf_line_parser import VCFLineSVPopulation
from upsetplot import plot
from matplotlib import pyplot
from upsetplot import from_memberships
from DataClasses import *
from dataclasses import dataclass
from functions_variables_single import *

class GenomeChartDataGenerator:
    def __init__(self, input_file_path, output_file_path):
        self.input_file_path = input_file_path
        self.output_file_path = output_file_path

    def allele_frequency_chart_generator(self):
        DEL_LIST = []
        INS_LIST = []
        INV_LIST = []
        DUP_LIST = []

        with open(self.input_file_path, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                obj = VCFLineSVPopulation(line)
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

        NUMBER_SAMPLES = len(DEL_LIST[0][:100])
        AF_single_sample_flag = int(NUMBER_SAMPLES == 1)

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
            x_label = 'AF_sample_' + str(i + 1) if not AF_single_sample_flag else 'AF'
            plt.xlabel(x_label)
            plt.ylabel('Count (log)')
            plt.title('Variant Frequency Spectrum')
            plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1), prop={'size': 5})

            plt.savefig(self.output_file_path + "AF_sample_" + str(i + 1), dpi=800)
            plt.close()



    
    def samples_sv_numbers(self):
        sum_GT=0
        f_sv_samples=open(self.output_file_path+"tmp.txt","w")
        with open(self.input_file_path, "r") as f:
            lines = f.readlines()
    
            for line in lines:
                if line.startswith("#"):
                    continue
                obj = VCFLineSVPopulation(line)
                if obj.FILTER=="PASS":
                    f_sv_samples.write(obj.SUPP_VEC+"\n")
    
        cmd = f"sort {self.output_file_path}tmp.txt | uniq -c > {self.output_file_path}sv_sample_results.txt"
        subprocess.run(cmd, shell=True)
        cmd= f"sed -i 's/^ *//' {self.output_file_path}sv_sample_results.txt"
        subprocess.run(cmd, shell=True)
        data = []
        data_set = []
        with open(f"{self.output_file_path}sv_sample_results.txt","r") as f:
                
            lines=f.readlines()
            for line in lines:
                elements = line.split()
                tmp_line = elements[-1].rstrip("\n")
                data_set.append(tmp_line)
                data.append(int(elements[0]))
            data_list=[]
            for item in data_set:
                tmp_list=[]
                if int(item) !=0:
                    for i in range(len(item)):
                        if item[i] != "0":
                            tmp="sample_"+str(i+1)
                            tmp_list.append(tmp)
                    data_list.append(tmp_list)
        
                elif int(item) ==0:
                    data_list.append([])
        
        example = from_memberships(data_list,data=data)
        upset = plot(example, facecolor="black", other_dots_color=.4,shading_color=.1)
    
        for patch in upset['intersections'].patches:
            upset['intersections'].annotate(text=patch.get_height(), xy=(patch.get_x() + patch.get_width() / 2, patch.get_height()), ha='center', va='bottom', rotation='vertical', fontsize=6, xytext=(0, +5), textcoords='offset points')
    
        for patch in upset['intersections'].patches:
            print(patch)
        upset['intersections'].bar_color = 'blue'
        upset['intersections'].bar_alpha = 0.7
        pyplot.savefig(f"{self.output_file_path}sample_upset.png",dpi=800, edgecolor="white")






