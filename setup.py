#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = ["biopython>=1.73", "python_codon_tables>=0.1.8", 'Bio']

setup_requirements = []

test_requirements = []

setup(
    author="Brian D. Weitzner",
    author_email="bweitzner@lyellbio.com",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    description="Amino acid reverse translation and DNA optimization tool based on species-specific codon-use distributions.",
    entry_points={
        "console_scripts": ["codon_harmony=codon_harmony.codon_harmony:main"]
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="codon_harmony",
    name="codon_harmony",
    packages=find_packages(exclude=["tests"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    python_requires=">=3.6",
    version="1.0.0",
    zip_safe=False,
)
