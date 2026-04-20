from setuptools import find_packages, setup


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
            "ldsc=ldsc.cli:main",
        ]
    },
    install_requires=[
        "bitarray>=2.5,<2.6",
        "scipy>=1.11,<1.12",
        "numpy>=1.24,<1.25",
        "pandas>=2.0,<2.1",
        "pyarrow>=16,<17",
        "pyliftover>=0.4,<0.5",
    ],
)
