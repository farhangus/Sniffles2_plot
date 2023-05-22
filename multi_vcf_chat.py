# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:17:54 2023

@author: HGSC
"""
import matplotlib.pyplot as plt
import numpy as np

from vcf_line_parser import VCFLineSVPopulation
f=open('hg003-GRCh37-ONT.vcf')
lines=f.readlines()
DEL_LIST=[]
INS_LIST=[]
INV_LIST=[]
for line in lines[3:]:
    if line[0] != "#":
        obj=VCFLineSVPopulation(line)
        tmp_sum=(int(obj.samples_DV[0])+int(obj.samples_DR[0]))
        if tmp_sum == 0:
            tmp_sum=0.01
#        print(tmp_sum)
#        print(obj.samples_DR)
#        print(obj.samples_DV)
#        print(obj.samples_DR)
#        print(obj.N_SUPP_VEC)
        allele_frq=(int(obj.samples_DV[0])/tmp_sum)
        if obj.SVTYPE=="DEL":
            DEL_LIST.append(allele_frq)
        if obj.SVTYPE=="INS":
            INS_LIST.append(allele_frq)
        if obj.SVTYPE=="INV":
            INV_LIST.append(allele_frq)
# Define the bin size and number of bins
bin_size = 0.04
num_bins = int(1 / bin_size)

# Create the histogram for each list
plt.hist([DEL_LIST, INS_LIST, INV_LIST], bins=num_bins, range=(0, 1), label=['DEL', 'INS', 'INV'], alpha=0.7, edgecolor='black')
plt.yscale('log')

# Set the x-axis and y-axis labels
plt.xlabel('AF')
plt.ylabel('Count (log)')
# Set the title of the plot
plt.title('Site Frequency Spectrum')
plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1),prop={'size': 5})
plt.savefig("test.jpg",dpi=800)

# Add a legend

