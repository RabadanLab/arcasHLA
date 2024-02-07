import io
import os
from setuptools import find_packages, setup

MAJOR = "0"
MINOR = "0.1"
VERSION = "{}.{}".format(MAJOR, MINOR)

with io.open("README.md", "r", encoding="utf-8") as f:
    README = f.read()

setup(
    name="arcas-hla",
    version=VERSION,
    description="High resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.",
    long_description=README,
    long_description_content_type="text/x-rst",
    namespace_packages=[],
    package_dir={"": "src"},
    packages=find_packages(where="src", exclude=["test*"]),
    include_package_data=True,
    python_requires=">=3.6,<4",
)