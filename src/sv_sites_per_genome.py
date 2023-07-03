# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import os
import pandas as pd
import numpy as np
import seaborn as sns
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
        self.data={}
        self.del_list = []
        self.dup_list = []
        self.ins_list = []
        self.inv_list = []
        self.bnd_list = []
        self.cnv_list = []
        self.other_list = []
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
                        if obj.SVTYPE == "DEL":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.del_list.append(list(integer_list))
                        elif obj.SVTYPE == "INS":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.ins_list.append(list(integer_list))
                        elif obj.SVTYPE == "DUP":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.dup_list.append(list(integer_list))
                        elif obj.SVTYPE == "INV":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.inv_list.append(list(integer_list))
                        elif obj.SVTYPE == "BND":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.bnd_list.append(list(integer_list))
                        elif obj.SVTYPE == "CNV":
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.cnv_list.append(list(integer_list))
                        else:
                            integer_list = [int(x) for x in obj.SUPP_VEC]
                            self.other_list.append(list(integer_list))
                            
                            


        df=pd.DataFrame(self.del_list)        
        column_sums = df.sum()
        series_list = column_sums.tolist()
        df.truncate(before=-1)
        self.data['del']=series_list
        
        df=pd.DataFrame(self.ins_list)        
        column_sums = df.sum()
        series_list = column_sums.tolist()
        df.truncate(before=-1)
        self.data['ins']=series_list
        
        # df=pd.DataFrame(self.inv_list)        
        # column_sums = df.sum()
        # series_list = column_sums.tolist()
        # df.truncate(before=-1)
        # self.data['inv']=series_list
        
        # df=pd.DataFrame(self.dup_list)        
        # column_sums = df.sum()
        # series_list = column_sums.tolist()
        # df.truncate(before=-1)
        # self.data['dup']=series_list
        
        # df=pd.DataFrame(self.bnd_list)        
        # column_sums = df.sum()
        # series_list = column_sums.tolist()
        # df.truncate(before=-1)
        # self.data['bnd']=series_list
        
        # df=pd.DataFrame(self.cnv_list)        
        # column_sums = df.sum()
        # series_list = column_sums.tolist()
        # df.truncate(before=-1)
        # self.data['cnv']=series_list
        
        # df=pd.DataFrame(self.other_list)        
        # column_sums = df.sum()
        # series_list = column_sums.tolist()
        # df.truncate(before=-1)
        # self.data['other']=series_list
        # print(self.data)
        
        df_final = pd.DataFrame(self.data)
        
        # Create the violin plot
        plt.figure(figsize=(10, 6))
        ax = sns.violinplot(data=df_final)
        
        plt.xlabel('SV', fontweight='bold', color='black')
        plt.ylabel('Count',  fontweight='bold', color='black')
        plt.title('SV sites per genome',fontweight='bold', color='black')
        
        # Annotate mean, best, and worst model values
        for i, dataset in enumerate(df_final.columns):
            mean_value = df_final[dataset].mean()
            max_value = df_final[dataset].max()
            min_value = df_final[dataset].min()
            
            ax.text(i, mean_value, f'mean:{mean_value:.0f}', ha='center', va='bottom', color='black')
            ax.text(i, max_value, f'max:{max_value:.0f}', ha='center', va='bottom', color='black')
            ax.text(i, min_value, f'min:{min_value:.0f}', ha='center', va='top', color='black')
        
        # Remove upper and right borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Make model names bold
        ax.set_xticklabels(ax.get_xticklabels(), fontweight='bold')
        
        
        plt.tight_layout()    
        plt.savefig(self.output_file('sv_site_per_genome.jpg'))

       #  df = pd.DataFrame(data)
       #  df = df.set_index('Model')
        
       #  # Create the violin plot
       #  plt.figure(figsize=(10, 6))
       #  ax = sns.violinplot(data=df)

       #  plt.savefig('kirekhar.jpg')

       #  exit()
       #  df = pd.DataFrame(self.del_list)
       #  df.appned(self.ins_list)
       #  column_sums = df.sum()
       #  series_list = column_sums.tolist()
       #  print(series_list)
       #  mean_value = np.mean(series_list)
       # # Create the violin plot
       #  sns.violinplot(data=series_list)

       #  plt.annotate('+', xy=(0.5, 0.5), xytext=(0.5, 0.5), xycoords='axes fraction',
       #       textcoords='axes fraction', fontsize=40, ha='center', va='center')
       #  plt.title("Violin Plot")
       #  plt.xlabel("Values")
       #  plt.ylabel("Density")
       #  plt.savefig('kirekhar.jpg')
