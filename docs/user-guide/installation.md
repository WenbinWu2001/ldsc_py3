# Installation

Last updated on: 2026-10-27

LDSC3 is installed from its Git repository into a conda environment, as an editable package so that source updates take effect immediately.
Choose the stable `main` branch for routine use, or the `restructure` branch for the current development version.

## Stable version

Clone the default `main` branch into a directory `ldsc3`:

```bash
git clone https://github.com/WenbinWu2001/ldsc_py3.git ldsc3
cd ldsc3
```

Create the conda environment and install the package with development extras:

```bash
conda env create -f environment.yml
conda activate ldsc3
python -m pip install -e ".[dev]"
```

Activate the environment before running LDSC3 commands or scripts in later sessions:

```bash
conda activate ldsc3
```

To update, pull the latest changes from `main`:

```bash
cd ldsc3
git pull origin main
conda activate ldsc3
```

If `environment.yml` or the package dependencies changed, also update the environment and reinstall the package:

```bash
conda env update -f environment.yml --prune
python -m pip install -e ".[dev]"
```

## Development version

Clone the `restructure` branch into a directory `ldsc3-dev`:

```bash
git clone --branch restructure https://github.com/WenbinWu2001/ldsc_py3.git ldsc3-dev
cd ldsc3-dev
```

Create a separate development environment and install the package with development extras:

```bash
conda env create -f environment.yml -n ldsc3-dev
conda activate ldsc3-dev
python -m pip install -e ".[dev]"
```

Activate it in later sessions with `conda activate ldsc3-dev`.
To update, run `git pull origin restructure` in `ldsc3-dev`, and if `environment.yml` or the dependencies changed, run:

```bash
conda env update -f environment.yml -n ldsc3-dev --prune
python -m pip install -e ".[dev]"
```

## Python versions and extras

LDSC3 supports Python 3.11 to 3.13.
The base installation includes NumPy, pandas, SciPy, PyArrow, Matplotlib, threadpoolctl and joblib; threadpoolctl limits native numerical libraries inside workers, and joblib supplies the detected CPU budget ([parallel execution](../reference/parallel-execution.md)).
Matplotlib is installed by default, but plotting never runs automatically and is imported only when a figure is made.
BED and gene-list interval projection uses NumPy from the base installation; no extra package or external BED tool is needed.

| Extra | Installs | Needed for |
| --- | --- | --- |
| `.[plink]` | `bitarray` | LD-score computation from PLINK files |
| `.[liftover]` | `pyliftover` | Chain-file liftover in `munge-sumstats` and `build-r2-panel` |
| `.[dev]` | the PLINK and liftover extras plus pytest | Development and the test suite |

## Check the installation

```bash
ldsc --help
python -m ldsc --help
```

Both commands list the same subcommands; continue with the [guided tutorial](guided-tutorial.md).
