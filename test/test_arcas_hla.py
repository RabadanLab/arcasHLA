# Test partial and whole HLA typing pipelines
import pytest
import subprocess
import json
from os.path import dirname, abspath


ROOT_DIR = dirname(dirname(abspath(__file__)))

@pytest.fixture(scope="session", autouse=True)
def set_reference_version():
    """
    Fetch IMGT/HLA database version 3.24.0 before test suite
    """
    reference_cmd = f"{ROOT_DIR}/arcasHLA reference --version 3.24.0"
    subprocess.run(reference_cmd.split())


@pytest.fixture(scope="session", autouse=True)
def extract_reads():
    """
    Extract reads before typing tests
    """
    extract_cmd = f"{ROOT_DIR}/arcasHLA extract test/test.bam -o test/output -t 8 -v"
    subprocess.run(extract_cmd.split())


def test_whole_allele_typing():
    whole_typing_cmd = (
        f"{ROOT_DIR}/arcasHLA genotype test/output/test.extracted.1.fq.gz "
        f"test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o test/output -t 8 -v"
    )
    subprocess.run(whole_typing_cmd.split())

    output_file = "./test/output/test.genotype.json"
    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*39:01:01", "B*07:02:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*06:09:01", "DQB1*02:02:01"],
        "DRB1": ["DRB1*10:01:01", "DRB1*14:02:01"]
    }
    with open(output_file, "r") as f:
        output = json.load(f)
    assert(output == expected_output)


def test_partial_allele_typing():
    partial_typing_cmd = (
        f"{ROOT_DIR}/arcasHLA partial test/output/test.extracted.1.fq.gz "
        f"test/output/test.extracted.2.fq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -G test/output/test.genotype.json "
        f"-o test/output -t 8 -v"
    )
    subprocess.run(partial_typing_cmd.split(), shell=True)

    output_file = "./test/output/test.partial_genotype.json"
    expected_output = {
        "A": ["A*01:01:01", "A*03:01:01"],
        "B": ["B*07:02:01", "B*39:39:01"],
        "C": ["C*08:01:01", "C*01:02:01"],
        "DPB1": ["DPB1*14:01:01", "DPB1*02:01:02"],
        "DQA1": ["DQA1*02:01:01", "DQA1*05:03"],
        "DQB1": ["DQB1*02:02:01", "DQB1*06:04:01"],
        "DRB1": ["DRB1*03:02:01", "DRB1*14:02:01"]
    }
    with open(output_file, "r") as f:
        output = json.load(f)
    assert(output == expected_output)
