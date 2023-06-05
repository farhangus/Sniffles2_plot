# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 10:48:59 2023

@author: HGSC
"""

import sys
import os

def get_help():
    help_message = '''
    Usage: ./single_visualizer_pipeline.py [OPTIONS]
    Options:
      -h  Display this help message.
      -i  Specify the path to the VCF files directory (accepts both single and multi samples VCF files).
    Example:
      ./pipeline.py -i <vcfs/>
    '''
    print(help_message)
    sys.exit(1)

if len(sys.argv) == 1:
    get_help()

input_file_path = None
i = 1

while i < len(sys.argv):
    arg = sys.argv[i]
    if arg == '-h':
        get_help()
    elif arg == '-i':
        if i + 1 < len(sys.argv):
            input_file_path = sys.argv[i + 1]
            i += 1
    elif arg == '-o':
        if i + 1 < len(sys.argv):
            output_file_path = sys.argv[i + 1]
            i += 1
        else:
            get_help()
    else:
        get_help()
    i += 1

print(input_file_path)
if os.path.isfile(input_file_path):
        file_name = os.path.splitext(input_file_path)[0]
        print(file_name)
        directory = os.path.dirname(os.path.abspath(input_file_path))
        os.makedirs(directory, exist_ok=True)
        print("Created directory:", directory)
        print(directory)
        os.system(f"python3 vcf_visulaizer.py -i {input_file_path} -o {output_file_path}/")
else:
    for file_name in os.listdir(input_file_path):
        file_path = os.path.join(input_file_path, file_name)
        if os.path.isfile(file_path):
            file_name = os.path.splitext(file_name)[0]
            directory_path = os.path.join(input_file_path, file_name)
            os.makedirs(directory_path, exist_ok=True)
            print("Created directory:", directory_path)
            os.system(f"python3 vcf_visulaizer.py -i {file_path} -o {directory_path}/")
