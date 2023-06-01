#!/bin/bash

get_help() {
    cat <<EOF
     "Usage: ./single_visualizer_pipeline.sh [OPTIONS]"
     "Options:"
     "  -h  Display this help message."
     "  -i  Specify the path to the VCF files directory."
     "Example:"
     "  ./single_visualizer_pipeline.sh -i <tmp_single_vcf/>"
    
EOF
    exit 1
}

if [ "$#" -eq 0 ]
then
    get_help
fi

while [ "$#" -gt 0 ]
do
    case "$1" in
        -[h])get_help;;
        -[i])INPUT_FILE_PATH="$2";shift;;
        *) get_help;;
    esac
    shift
done


# Iterate through the files in the directory
for file_path in "$INPUT_FILE_PATH"/*; do
    if [[ -f "$file_path" ]]; then
        # Extract the file name without extension
        file_name=$(basename "$file_path")
        file_name="${file_name%.*}"
            # Create the directory
            rm -f -r "$INPUT_FILE_PATH/$file_name"
            mkdir "$INPUT_FILE_PATH/$file_name"
            echo "Created directory: $INPUT_FILE_PATH/$file_name"
            python3 vcf_visulaizer.py -i "$file_path" -o "$INPUT_FILE_PATH/$file_name/"
    fi
done
