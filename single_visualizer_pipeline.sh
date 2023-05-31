#!/bin/bash

# Specify the directory path
directory_path="$1"

# Iterate through the files in the directory
for file_path in "$directory_path"/*; do
    if [[ -f "$file_path" ]]; then
        # Extract the file name without extension
        file_name=$(basename "$file_path")
        file_name="${file_name%.*}"
            # Create the directory
            rm -f -r "$directory_path/$file_name"
            mkdir "$directory_path/$file_name"
            echo "Created directory: $directory_path/$file_name"
            python3 vcf_visulaizer.py -i "$file_path" -o "$directory_path/$file_name/"
    fi
done
