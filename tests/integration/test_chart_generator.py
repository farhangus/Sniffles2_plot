from sniffles2_plot.cli import multi_visualizer
from pathlib import Path
from os import getcwd, listdir

def test_generate_should_output_expected_files(tmp_path:Path) -> None:
    #arrange

    input_vcf_file = Path(__file__).parent.parent / "fixtures/hg002-trio_merge.vcf" 
    output_path = tmp_path/"generated"
    output_path.mkdir()
    # act
    multi_visualizer(input_vcf_file,output_path)
    actual_output_file_names=listdir(output_path)
    #assert
    assert len(actual_output_file_names) == 8
