# -*- coding: utf-8 -*-
"""
Created on Wed May 17 10:40:52 2023

@author: HGSC
"""
import matplotlib.pyplot as plt
import numpy as np
from vcf_line_parser import VCFLineSV

ranges = [(50, 100), (100,200),(200,300),(300,400),(400,600),(600,800),(800,1000),(1000,2500),(2500,5000),(5000,10000000)]
x_labels = ['50-100', '100-200', '200-300','300-400','400-600','600-800','800-1k','1k-2.5k','2.5k-5k','>5k']

def count_numbers_in_ranges(numbers, ranges):
    counts = [0] * len(ranges)  # Initialize a count for each range
    
    for number in numbers:
        for i, range_ in enumerate(ranges):
            if range_[0] <= number <= range_[1]:
                counts[i] += 1  # Increment the count for the corresponding range
                break  # Exit the inner loop once the number is found in a range
    
    return counts


def sv_size_type_chart(l1,l2,l3,l4,l5,output_chart):
    
    fig, ax = plt.subplots()
    x = np.arange(len(l1))
    # Width of each bar
    bar_width = 0.14
   # Creating the bar plots
    plt.bar(x, l1, width=bar_width, label='DEL')
    plt.bar(x + bar_width, l2, width=bar_width, label='INS')
    plt.bar(x + 2*bar_width, l3, width=bar_width, label='DUP')
    plt.bar(x + 3*bar_width, l4, width=bar_width, label='INV')
    plt.bar(x + 4*bar_width, l5, width=bar_width, label='BND')

    # Adding labels and title
    plt.xlabel('SizeBin')
    plt.ylabel('Count')
    plt.title('SV Size/Type Distribution')
    plt.xticks(x + bar_width/5, x)
    ax.set_xticklabels(x_labels,rotation=45)
    plt.tight_layout()  # Adjust the padding to accommodate all labels
    plt.legend(title='svtype')
    # Displaying the plot
    plt.savefig(output_chart,dpi=800)
    
def lenght_var_count_chart(output_name,del_flag,v_s_list,min_val,max_value,chart_pos,ytick_flag,xtick_locs_list,plt_titl,subplt_title):

    plt.subplot(1,6, chart_pos)
    plt.subplots_adjust(wspace=0.05, hspace=0.05,bottom=0.3)

    minimum = min_val
    maximum = max_value
    selected_values = [x for x in v_s_list if minimum <= x <= maximum]
    # Determine the bin size
    bin_size = round((max(selected_values) - min(selected_values))/1000)
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
    f=open(input_vcf_file)
    lines=f.readlines()
    DEL=[]
    INS=[]
    DUP=[]
    INV=[]
    BND=[]
    DEL_NEGATIVE=[]
    for line in lines:
        if line[0] != "#":
            obj=VCFLineSV(line)
            if obj.SVTYPE=="DEL":
                DEL.append(abs(obj.SVLEN))
            elif obj.SVTYPE=="INS":
                INS.append(obj.SVLEN)
            elif obj.SVTYPE=="INV":
                INS.append(obj.SVLEN)
            elif obj.SVTYPE=="DUP":
                INS.append(obj.SVLEN)
            elif obj.SVTYPE=="BND":
                INS.append(obj.SVLEN)
    return DEL,INS,DUP,INV,BND




