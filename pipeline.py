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
      -i  Specify the path to the VCF files directory.
    Example:
      ./single_visualizer_pipeline.py -i <tmp_single_vcf/>
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
        else:
            get_help()
    else:
        get_help()
    i += 1

if input_file_path is None:
    get_help()

for file_name in os.listdir(input_file_path):
    file_path = os.path.join(input_file_path, file_name)
    if os.path.isfile(file_path):
        file_name = os.path.splitext(file_name)[0]
        directory_path = os.path.join(input_file_path, file_name)
        os.makedirs(directory_path, exist_ok=True)
        print("Created directory:", directory_path)
        os.system(f"python3 vcf_visulaizer.py -i {file_path} -o {directory_path}/")
