import argparse
import os
from typing import Tuple
from sniffles2_plot.cli.generate_charts import generate_charts 

def _parse_arguments() -> Tuple[str, str]:
    parser = argparse.ArgumentParser(description="Snfiffle2 plot generator")
    parser.add_argument(
        "-i",
        "--input",
        help="Specify the path to the VCF files directory (accepts both single and multi samples VCF files)",
        required=True,
    )
    parser.add_argument(
        "-o", "--output", help="Specify the path to the VCF output file"
    )
    args = parser.parse_args()
    return args.input, args.output


def _ensure_output_directory_exists(path:str) -> None:
        if os.path.exists(path) and not os.path.isdir(path):
            raise IOError("The given path is not a directory.")
        if not os.path.exists(path):
            os.mkdir(path)


def entry_point():
    input_file_path, output_directory_path = _parse_arguments()

    if os.path.isfile(input_file_path):
        _ensure_output_directory_exists(output_directory_path)
        generate_charts(input_file_path, output_directory_path)
    else:
        for file_name in os.scandir(input_file_path):
            if file_name.name.lower().endswith(".vcf"):
                directory_path = os.path.join(input_file_path, os.path.splitext(file_name.name)[0])
                os.makedirs(directory_path, exist_ok=True)
                generate_charts(file_name, directory_path)
                