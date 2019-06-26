#! /usr/bin/env python

import os

from distutils.core import setup
from setuptools.command.test import test as TestCommand

NAME = 'simcalc'


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

install_requires = read('requirements.txt')

extras_require = {}
# Notebook dependencies for plotting
extras_require['plotting'] = read('requirements_plotting.txt')

setup(
    name=NAME,
    install_requires = install_requires,
    extras_require = extras_require,
    version="1.0",
    author="Sander Keemink",
    author_email="swkeemink@scimail.eu",
    description="Simulate calcium data",
    url="https://github.com/rochefort-lab/SimCalc",
    download_url="NA",
    package_dir={NAME: "./simcalc"},
    packages=[NAME],
    license="GNU",
    long_description=read('README.rst'),
    classifiers=[
        "License :: GNU",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering"
    ],
)
