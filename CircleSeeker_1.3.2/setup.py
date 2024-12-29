"""
CircleSeeker: Advanced Circular DNA Analysis Tool.

A comprehensive bioinformatics package designed for detecting and analyzing
circular DNA from PacBio HiFi sequencing data. This tool provides a complete
pipeline from raw sequence processing to detailed result visualization.

Features:
- Specialized algorithms for circular DNA detection
- Multi-threaded processing for high performance
- Comprehensive reporting and visualization
- Support for various circular DNA types

For more information, visit: https://github.com/leoxqy/CircleSeeker

Copyright (c) 2024 CircleSeeker Team
"""

from setuptools import setup, find_packages

setup(
    name='circle_seeker',
    version='1.3.2',
    packages=find_packages(),
    description='Circle-Seeker: A specialized bioinformatics tool designed for advanced circular DNA analysis using PacBio HiFi sequencing data.',
    author='Yaoxin Zhang & Leo Xinqi Yu',
    author_email='yxzhang@ncgr.com, leoxqy@hotmail.com',
    url='https://github.com/leoxqy/CircleSeeker',
    license='GPL-2.0',
    install_requires=[
        'pandas',
        'numpy',
        'pysam',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'CircleSeeker=circle_seeker.CircleSeeker:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: OS Independent',
    ],
    include_package_data=False,
    python_requires='>=3.6',
)
