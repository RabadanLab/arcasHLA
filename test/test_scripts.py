"""
Test the individual scripts for expected output.
"""

import json
import subprocess


def test_whole_allele_typing(repo_root, extract_reads, tmp_path):
    output_dir = str(tmp_path)

    whole_typing_cmd = (
        f"{repo_root}/arcasHLA genotype {extract_reads[0]} "
        f"{extract_reads[1]} -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o {output_dir} -t 8 -v"
    )
    subprocess.run(whole_typing_cmd.split(), check=True)

    output_file = f"{output_dir}/test.genotype.json"

    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*39:01:01", "B*07:02:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*02:02:01", "DQB1*06:09:01"],
        "DRB1": ["DRB1*10:01:01", "DRB1*14:02:01"],
    }
    with open(output_file, "r") as f:
        output = json.load(f)

    # Convert lists in output and expected output to sets since order does not matter for test
    for key in output:
        output[key] = set(output[key])
        expected_output[key] = set(expected_output[key])
    assert output == expected_output


def test_partial_allele_typing(repo_root, extract_reads, tmp_path):
    output_dir = str(tmp_path)

    partial_typing_cmd = (
        f"{repo_root}/arcasHLA partial {extract_reads[0]} "
        f"{extract_reads[1]} "
        "-g A,B,C,DPB1,DQB1,DQA1,DRB1 "
        f"-G {repo_root}/test/expected_output/test.genotype.json "
        f"-o {output_dir} "
        "-t 8 "
        "-v"
    )
    subprocess.run(partial_typing_cmd.split(), check=True)

    output_file = f"{output_dir}/test.partial_genotype.json"
    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*07:02:01", "B*39:39:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*06:04:01", "DQB1*02:02:01"],
        "DRB1": ["DRB1*03:02:01", "DRB1*14:02:01"],
    }
    with open(output_file, "r") as f:
        output = json.load(f)

    # Convert lists in output and expected output to sets since order does not matter for test
    for key in output:
        output[key] = set(output[key])
        expected_output[key] = set(expected_output[key])
    assert output == expected_output


def test_reference_customization(repo_root, tmp_path):
    genotype_result_path = f"{repo_root}/test/expected_output/test.genotype.json"
    custom_ref_out_path = f"{str(tmp_path)}/custom_reference"

    custom_ref_command = [
        f"{repo_root}/arcasHLA",
        "customize",
        "--genotype",
        genotype_result_path,
        "-o",
        custom_ref_out_path,
    ]

    customize_result = subprocess.run(custom_ref_command, capture_output=True)

    assert customize_result.returncode == 0
