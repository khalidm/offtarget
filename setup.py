#!/usr/bin/env python

from setuptools import setup

setup(
    name='offtarget',
    version='0.0.1',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['offtarget'],
    package_data={},
    entry_points={
        'console_scripts': ['offtarget = offtarget.offtarget:main']
    },
    url='https://github.com/bjpop/offtarget',
    license='LICENSE.txt',
    description='Plot offtarget reads for targeted sequence data',
    #long_description=open('README.txt').read(),
    install_requires=[
        "pysam >= 0.8.1",
        "matplotlib",
        "datrie",
        "biopython"
    ],
    classifiers=[
    ],
)
