#!/usr/bin/env python
from MUTANT import __version__
from setuptools import setup, find_packages

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="MUTANT",
    version=version,
    long_description=__doc__,
    url="https://github.com/Clinical-Genomics/MUTANT",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': ['MUTANT=MUTANT.cli:root'],
    },
)