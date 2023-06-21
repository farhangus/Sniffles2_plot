# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
from vcf_line_parser import VCFLineSVPopulation



class VariantCount:
    """generate plots for variant count all SVs"""
    def __init__(self, input_file_path, output_directory):
        self.input_file_path = input_file_path
        self.output_directory = output_directory

    def output_file(self, filename):
        """returnthe output file name and filepath"""
        return os.path.join(self.output_directory,filename)

    def variant_count_chart_generator(self):
        """generate the allele frequency plots"""
        self.del_count=0
        self.dup_count=0
        self.ins_count=0
        self.inv_count=0
        self.bnd_count=0
        self.cnv_count=0
        self.other_count=0
        test_sv=[]
        with open(self.input_file_path, "r",encoding="utf-8") as f:
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
                    self.dup_count += 1
                elif sv_type == "INV":
                    self.ins_count += 1
                elif sv_type == "DUP":
                    self.dup_count += 1
                elif sv_type == "BND":
                    self.bnd_count += 1 
                elif sv_type == "CNV":
                    self.cnv_count += 1 
                else:
                    self.other_count +=1
        labels=["DEL","DUP","INS","INV","BND","CNV","OTHER"]
        frequencies = [self.del_count, self.dup_count, self.ins_count, 
                       self.inv_count, self.dup_count, self. bnd_count,
                       self.cnv_count, self.other_count]
        # Create the bar plot
        print(frequencies)
        plt.bar(frequencies,frequencies)
        
        # Add labels and title
        plt.xlabel("Types")
        plt.ylabel("Frequencies")
        plt.title("Frequency of Variant Types")
        
        # Show the plot
        plt.savefig(self.output_file('variant_count.jpg'))
                


    




