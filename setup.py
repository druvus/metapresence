from setuptools import setup, find_packages

setup(
    name="metapresence",
    version="1.0.0",
    description="Metrics for evaluating the distribution of aligned reads onto nucleotidic sequences",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pysam",
    ],
    extras_require={
        "plotting": ["matplotlib"],
    },
    entry_points={
        "console_scripts": [
            "metapresence=metapresence.cli:main",
            "metapresence-merge=metapresence.cli:merge",
            "metapresence-parse=metapresence.cli:parse",
        ],
    },
    python_requires=">=3.6",
)