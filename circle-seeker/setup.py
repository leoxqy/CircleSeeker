from setuptools import setup, find_packages

setup(
    name='circle-seeker',
    version='1.0.0',
    packages=find_packages(),
    description='Circle Seeker: a package for advanced circular DNA analysis.',
    author='Your Name',
    author_email='your.email@example.com',
    url='<your-repository-url>',
    install_requires=[
        'minimap2',
        'samtools',
        'blast',
        'bedtools',
        'seqkit',
        'tidehunter',
        'pandas',
        'numpy',
        'pysam',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'CircleSeeker=circle-seeker.CircleSeeker:main',
        ],
    }
)
