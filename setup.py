from setuptools import find_packages, setup


INSTALL_REQUIRES = [
    "scipy>=1.15,<1.18",
    "numpy>=2.1,<2.5",
    "pandas>=2.3,<3.0",
    # Required for canonical parquet I/O, row-group pruning, and LDSC schema metadata.
    "pyarrow>=21,<24",
]

PLINK_REQUIRES = [
    "bitarray>=3.0,<4",
]

BED_REQUIRES = [
    "pybedtools>=0.12,<0.13",
]

LIFTOVER_REQUIRES = [
    "pyliftover>=0.4.1,<0.5",
]

TEST_REQUIRES = [
    "pytest>=8,<10",
]

ALL_REQUIRES = PLINK_REQUIRES + BED_REQUIRES + LIFTOVER_REQUIRES


setup(
    name="ldsc",
    version="2.0b",
    description="LD Score Regression (LDSC)",
    url="http://github.com/abrantesas/ldsc_py3",
    author="Anthony Abrantes, Brendan Bulik-Sullivan and Hilary Finucane",
    author_email="antshaabr@gmail.com",
    license="GPLv3",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
    package_data={"ldsc": ["data/*.tsv.gz"]},
    entry_points={
        "console_scripts": [
            "ldsc=ldsc.cli:run_cli",
        ]
    },
    python_requires=">=3.11,<3.14",
    install_requires=INSTALL_REQUIRES,
    extras_require={
        "all": ALL_REQUIRES,
        "bed": BED_REQUIRES,
        "dev": ALL_REQUIRES + TEST_REQUIRES,
        "liftover": LIFTOVER_REQUIRES,
        "plink": PLINK_REQUIRES,
        "test": TEST_REQUIRES,
    },
)
