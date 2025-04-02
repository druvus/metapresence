from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="metapresence",
    version="1.0.0",
    description="Metrics for evaluating the distribution of aligned reads onto nucleotidic sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Original Authors",
    author_email="original_authors@example.com",
    url="https://github.com/original/metapresence",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pysam",
    ],
    extras_require={
        "plotting": ["matplotlib"],
        "interactive": ["plotly>=5.0.0", "scipy"],
        "all": ["matplotlib", "plotly>=5.0.0", "scipy"]
    },
    entry_points={
        "console_scripts": [
            "metapresence=metapresence.cli:main",
            "metapresence-merge=metapresence.cli:merge",
            "metapresence-parse=metapresence.cli:parse",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
)