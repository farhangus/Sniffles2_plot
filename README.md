# Sniffle plot

The Sniffle_plot is a Python program that generate plots for both single and mlulti samples VCF files.
Type of plots for single smaple VCF file:
    1- Varriant Frequency spectrum
    2- Genotype ferquency
    3- SV size &Type distribution
    4- Comparition of lenght of variants

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
     python3 sniffles2_plots.py -i <target directory>
     
For running the program for single vcf file:
     python3 sniffles2_plots.py -i <file_name> -o <output_folder>

