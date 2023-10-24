from os import getcwd, listdir
from pathlib import Path

from sniffles2_plot.cli import generate_multi_vcf_charts


def test_generate_should_output_expected_files(tmp_path: Path) -> None:
    # arrange

    input_vcf_file = Path(__file__).parent.parent / "fixtures/hg002-trio_merge.vcf"
    output_path = tmp_path / "generated"
    output_path.mkdir()
    # act
    generate_multi_vcf_charts(input_vcf_file, output_path)
    actual_output_file_names = listdir(output_path)
    # assert
    breakpoint()
    assert len(actual_output_file_names) == 7
