#
# setup.py
#
# Installation script for the LignoForge package.
#
import os
import sys
import setuptools

# Import the lengthy rich-text README as the package's long
# description:
root_dir = os.path.dirname(__file__)

with open(os.path.join(root_dir, "README.md"), "r") as fh:
	long_description = fh.read()


setuptools.setup(
    name="lignoforge",
    version="0.2.0",
    author="DAIMON Team",
    author_email="daimoners@gmail.com",
    description="Top-down lignin model generation from experimental input to atomistic polymers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/daimoners/lignoforge",
    project_urls={
        "Documentation": "https://lignoforge.readthedocs.io/en/latest/",
        "Source": "https://github.com/daimoners/lignoforge",
    },
    packages=setuptools.find_packages(),
    package_data={'': ['*.xlsx', '*.json']},
    include_package_data=True,
    python_requires=">=3.7",
    install_requires=[
        "matplotlib>=3.1.1",
        "numpy>=1.19.2",
        "scipy>=1.3.1",
        "pandas>=0.25.1",
        "openpyxl>=3.0.7",
        "pytest>=6.2.3",
        "jsonschema>=4.0.0",
        "networkx>=2.5",
        "pysmiles>=1.0.1",
        "rdkit-pypi>=2021.9.2.1"],
    entry_points={
        "console_scripts": [
            "lignoforge-chain = lignoforge.cli.build_chain:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Chemistry",
    ],
)