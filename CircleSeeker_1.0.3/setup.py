from setuptools import setup, find_packages

setup(
    name='circle_seeker',
    version='1.0.3',
    packages=find_packages(),
    description='Circle-Seeker: A specialized bioinformatics tool designed for advanced circular DNA analysis using PacBio HiFi sequencing data.',
    author='Yaoxin Zhang & Leo Xinqi Yu',
    author_email='yxzhang@ncgr.com, leoxqy@hotmail.com',
    url='https://github.com/leoxqy/CircleSeeker',
    license='GPL-2.0',
    install_requires=[
        'pandas',
        'numpy>=1.22',
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

