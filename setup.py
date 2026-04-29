from setuptools import find_packages, setup


INSTALL_REQUIRES = [
    "bitarray>=3.0,<4",
    "scipy>=1.15,<1.18",
    "numpy>=2.1,<2.5",
    "pandas>=2.3,<3.0",
    "pyarrow>=21,<24",
    "pybedtools>=0.12,<0.13",
    "pyliftover>=0.4.1,<0.5",
]

TEST_REQUIRES = [
    "pytest>=8,<10",
]


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
        "test": TEST_REQUIRES,
    },
)
