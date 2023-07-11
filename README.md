# Sniffle Plot

The Sniffle Plot is a Python program that generates plots for both single and multi-sample VCF files.

Type of plots for single sample VCF file:

1. Variant Frequency Spectrum
2. Genotype Frequency
3. SV Size & Type Distribution
4. Comparison of Length of Variants

Type of plots for multi smaples VCF file:
    1-Varriant Frequency spectrum
    2-Upset plot for sample intersection




## Requirements

The following dependencies are required to use the Genome Chart Data Generator:

- Python 3.x
- matplotlib
- upsetplot

## usage
For running the program for numbers of vcf files in a sepcific directory:
     python3 sniffles2_plots.py -i <VCF_files_folder>
     
For running the program for single vcf file:
     python3 sniffles2_plots.py -i <file_name> -o <output_folder>

