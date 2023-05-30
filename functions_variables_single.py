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
from DataClasses import VcfVariables
from dataclasses import dataclass

ranges = [(50, 100), (100,200),(200,300),(300,400),(400,600),(600,800),(800,1000),(1000,2500),(2500,5000),(5000,10000000)]
x_labels = ['50-100', '100-200', '200-300','300-400','400-600','600-800','800-1k','1k-2.5k','2.5k-5k','>5k']

def separate_lists(lst,num):
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

@dataclass
class GenomeChartData():
    points: list
    legend: str
    def count(self,labels):
        return [
            self.points.count(lbl)
            for lbl in labels
        ]
def genome_bar_chart(output_file_path,labels,*bars:GenomeChartData):
    print(output_file_path, labels, bars[0])
    x = np.arange(len(labels))
    # Define the width of each bar
    width = 0.2
    # Create the grouped bar plot
    plt.bar(x - width, bars[0].count(labels), width=width, label=bars[0].legend)
    plt.bar(x, bars[1].count(labels), width=width, label=bars[1].legend)
    # Add labels and title
    plt.xlabel("Genotype",size=6)
    plt.ylabel("Count",size=6)
    plt.title("Genotype_Frequency")
    # Set the x-axis tick labels
    plt.xticks(x, labels)
    # Add a legend
    plt.legend()
    plt.savefig(output_file_path,dpi=800)

def sv_size_type_chart(l1,l2,label_1,label_2,path_2_chart):
    fig, ax = plt.subplots()
    x = np.arange(len(l1))
    # Width of each bar
    bar_width = 0.40
   # Creating the bar plots
    plt.bar(x, l1, width=bar_width, label=label_1)
    plt.bar(x + bar_width, l2, width=bar_width, label=label_2)
#    plt.bar(x + 2*bar_width, l3, width=bar_width, label='DUP')
#    plt.bar(x + 3*bar_width, l4, width=bar_width, label='INV')
#    plt.bar(x + 4*bar_width, l5, width=bar_width, label='BND')

    # Adding labels and title
    plt.xlabel('SizeBin')
    plt.ylabel('Count')
    plt.title('SV Size/Type Distribution')
    plt.xticks(x + bar_width/5, x)
    ax.set_xticklabels(x_labels,rotation=45)
    plt.tight_layout()  # Adjust the padding to accommodate all labels
    plt.legend(title='svtype')
    # Displaying the plot
    plt.savefig(path_2_chart,dpi=800)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def lenght_var_count_chart(output_name,del_flag,v_s_list,min_val,max_value,chart_pos,ytick_flag,xtick_locs_list,plt_titl,subplt_title):

    plt.subplot(1,6, chart_pos)
    plt.subplots_adjust(wspace=0.05, hspace=0.05,bottom=0.3)

    minimum = min_val
    maximum = max_value
    selected_values = [x for x in v_s_list if minimum <= x <= maximum]
    # Determine the bin size
    bin_size = round((max(selected_values) - min(selected_values))/100)
    # Create bins based on the bin size
    bins = np.arange(min(selected_values), max(selected_values) + bin_size, bin_size)
    # Count the number of values that fall into each bin
    bin_counts, _ = np.histogram(selected_values, bins)
    # Display the bin plot

    plt.bar(bins[:-1], bin_counts, width=bin_size, align='edge')

    xtick_labels = [str(x) for x in xtick_locs_list]

    plt.xticks(xtick_locs_list, xtick_labels,rotation=90)
    if chart_pos==1:
        plt.ylabel('Counts')
    if chart_pos==3:
        plt.xlabel('\n\n        Lenght of Variants')

    plt.yscale('log')
    if ytick_flag:
        plt.yticks([1,10, 100,1000,10000])
    else:
        plt.yticks([1,10,100,1000,10000])
        plt.yticks([])
    plt.title(subplt_title,size=5)
    if del_flag:
        plt.gca().invert_xaxis()

    plt.savefig(output_name,dpi=1000)




def vcf_number_variants(input_vcf_file):
    with open(input_vcf_file,"r") as f:
        lines=f.readlines()
        vcf_variables=VcfVariables.new()

        for line in lines:
            if line[0] != "#":
                obj=VCFLineSV(line)
                if obj.SVTYPE=="DEL":
                    vcf_variables.DEL.append(abs(obj.SVLEN))
                    vcf_variables.DEL_GENOTYPE.append(obj.GENOTYPE)
                    vcf_variables.DEL_AF.append(obj.AF)
                elif obj.SVTYPE=="INS":
                    vcf_variables.INS.append(obj.SVLEN)
                    vcf_variables.INS_GENOTYPE.append(obj.GENOTYPE)
                    vcf_variables.INS_AF.append(obj.AF)
                elif obj.SVTYPE=="INV":
                    vcf_variables.INV.append(obj.SVLEN)
                    vcf_variables.INV_GENOTYPE.append(obj.GENOTYPE)
                    vcf_variables.INV_AF.append(obj.AF)
                elif obj.SVTYPE=="DUP":
                    vcf_variables.DUP.append(obj.SVLEN)
                    vcf_variables.DUP_GENOTYPE.append(obj.GENOTYPE)
                    vcf_variables.DEL_AF.append(obj.AF)
                elif obj.SVTYPE=="BND":
                    vcf_variables.BND.append(obj.SVLEN)
                    vcf_variables.BND_GENOTYPE.append(obj.GENOTYPE)
                    vcf_variables.BND_AF.append(obj.AF)
    return vcf_variables

def allele_frequency_chart_genrator(input_file_path,output_file_path):

    DEL_LIST = []
    INS_LIST = []
    INV_LIST = []
    DUP_LIST=[]
    with open(input_file_path, "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith("#"):
                continue
            obj = VCFLineSVPopulation(line)
            if obj.SVTYPE == "DEL":
                DEL_LIST.append(obj.samples_AF)
            elif obj.SVTYPE == "INS":
                INS_LIST.append(obj.samples_AF)
            elif obj.SVTYPE == "INV":
                INV_LIST.append(obj.samples_AF)
            elif obj.SVTYPE == "DUP":
                DUP_LIST.append(obj.samples_AF)

    
    NUMBER_SAMPLES=len(DEL_LIST[0][:100])
    AF_single_sample_flag = int(NUMBER_SAMPLES == 1)

    for i in range(NUMBER_SAMPLES):
        tmp_list_name_del="del_AF_"+str(i)
        tmp_list_name_del=[]
        tmp_list_name_del=separate_lists(DEL_LIST,i)
        tmp_list_name_ins="ins_AF_"+str(i)
        tmp_list_name_ins=[]
        tmp_list_name_ins=separate_lists(INS_LIST,i)
        tmp_list_name_inv="inv_AF_"+str(i)
        tmp_list_name_inv=[]
        tmp_list_name_inv=separate_lists(INV_LIST,i)
        tmp_list_name_dup="dup_AF_"+str(i)
        tmp_list_name_dup=[]
        tmp_list_name_dup=separate_lists(DUP_LIST,i)
        bin_size = 0.04
        num_bins = int(1 / bin_size)

        plt.hist([tmp_list_name_del, tmp_list_name_ins, tmp_list_name_inv,tmp_list_name_dup], bins=num_bins, range=(0, 1), label=['DEL', 'INS', 'INV','DUP'],
                  alpha=0.7, edgecolor='black')

        plt.yscale('log')
        x_label = 'AF_sample_' + str(i+1) if not AF_single_sample_flag else 'AF'
        plt.xlabel(x_label)
        plt.ylabel('Count (log)')
        plt.title('Varriant Frequency Spectrum')
        plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1), prop={'size': 5})

        plt.savefig(output_file_path+"AF_sample_"+str(i+1), dpi=800)
        plt.close()


def samples_sv_numbers(input_file_path,output_file_path):
    # DEL_LIST = []
    # INS_LIST = []
    # INV_LIST = []
    f_sv_samples=open(output_file_path+"tmp.txt","w")
    with open(input_file_path, "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith("#"):
                continue
            obj = VCFLineSVPopulation(line)
            f_sv_samples.write(obj.SUPP_VEC+"\n")
            # if obj.SVTYPE == "DEL":
            #     DEL_LIST.append(obj.samples_AF)
            # elif obj.SVTYPE == "INS":
            #     INS_LIST.append(obj.samples_AF)
            # elif obj.SVTYPE == "INV":
            #     INV_LIST.append(obj.samples_AF)

    cmd = f"sort {output_file_path}tmp.txt | uniq -c > {output_file_path}sv_sample_results.txt"
    subprocess.run(cmd, shell=True)
    cmd= f"sed -i 's/^ *//' {output_file_path}sv_sample_results.txt"
    subprocess.run(cmd, shell=True)
    data = []
    data_set = []
    f=open(f"{output_file_path}sv_sample_results.txt","r")
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



    # upset=plot(example,  facecolor="red", other_dots_color=.4,shading_color=.3)
    # for patch in upset['intersections'].patches:
    #     upset['intersections'].annotate(text=patch.get_height(), xy=(patch.get_x() + patch.get_width() /2, patch.get_height()), ha='center', va='bottom', rotation='vertical',fontsize=6,xytext=(0, +5), textcoords='offset points')
    # for index, row in enumerate(upset['totals']):
    #     upset['totals'].annotate(text=row.get_text(), xy=(index, 0), ha='center', va='bottom', fontsize=8)

    # plt.title('Intercection of samples',size=10)

    upset = plot(example, facecolor="black", other_dots_color=.4,shading_color=.1)

    # Rotate the counts vertically, adjust font size, and distance to the bar
    for patch in upset['intersections'].patches:
        upset['intersections'].annotate(text=patch.get_height(), xy=(patch.get_x() + patch.get_width() / 2, patch.get_height()), ha='center', va='bottom', rotation='vertical', fontsize=6, xytext=(0, +5), textcoords='offset points')

    # Show sum of the bottom chart in front of each row
    for patch in upset['intersections'].patches:
        print(patch)

    # Set the counts format

    # Customize the appearance of the bars
    upset['intersections'].bar_color = 'blue'
    upset['intersections'].bar_alpha = 0.7


    pyplot.savefig(f"{output_file_path}sample_upset.png",dpi=800, edgecolor="black")






