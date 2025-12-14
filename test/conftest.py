"""
Common fixtures, etc.
"""

import os
import pytest
import subprocess

from pathlib import Path


@pytest.fixture(scope="session")
def repo_root() -> str:
    """
    Provide the path to the repository root as defined as the parent to the dir in which
    this file lives.
    """
    return str(Path(__file__).parent.parent)


@pytest.fixture(scope="session", autouse=True)
def set_reference_version(repo_root):
    """
    Fetch IMGT/HLA database version 3.24.0 before test suite
    """
    reference_cmd = f"{repo_root}/arcasHLA reference --version 3.24.0"
    subprocess.run(reference_cmd.split(), check=True)


@pytest.fixture(scope="session")
def extract_reads(repo_root, tmp_path_factory):
    """
    Extract reads before typing tests
    """
    output_dir = str(tmp_path_factory.mktemp("extracted_reads"))
    extract_cmd = f"{repo_root}/arcasHLA extract test/test.bam -o {output_dir} -t 8 -v"
    subprocess.run(extract_cmd.split(), check=True)

    # Provide the individual extracted reads files.
    return [
        os.path.join(output_dir, outfile)
        for outfile in sorted(os.listdir(output_dir))
        if outfile.endswith(".fq.gz") and outfile.startswith("test")
    ]
