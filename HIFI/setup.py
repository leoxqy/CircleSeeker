from setuptools import setup, find_packages

setup(
    name="eccDNA_pipeline",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'biopython',
        'pandas',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'eccDNA_pipeline=eccDNA_pipeline.main:main',
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="A pipeline for eccDNA analysis",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/eccDNA_pipeline",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
