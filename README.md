[![DOI](https://zenodo.org/badge/654776766.svg)](https://zenodo.org/doi/10.5281/zenodo.13638032)

# GBM-Perturb

This repository contains analysis in Python and R for the output of _in vivo_ perturb-seq
experiments as described in [Liu et al, bioRxiv 2023](https://www.biorxiv.org/content/10.1101/2023.09.01.555831v3.full).

## Getting set up

### Installing dependencies
First, you'll want to install dependencies. Python and R dependencies are managed separately.

#### R
Dependencies are managed as an R project through [`renv`](https://rstudio.github.io/renv/articles/renv.html). Activate the project by installing `renv` and running:
```
renv::init()
renv::restore()
```
in the top-level project directory. This will install dependencies as specified in `renv.lock`.

#### Python
This repository uses Poetry for Python dependency management. [Install poetry](https://python-poetry.org/docs/), navigate to the
project directory, and run `poetry install`. This will install all packages according to the versions specified in `poetry.lock`.

### Directories
Next, you'll want to create a top-level `output` directory and a top-level `data` directory. The code in this repository assumes that
these directories (and sometimes subdirectories) exist. `output` generally contains output plots and data, while `data` contains
input metadata and raw data. We have chosen to leave most of our existing directory structure following these `data` and `output`
directories intact as examples; there's user flexibility on where things go.

## Navigating the repository
Here is the directory structure as used by us as developers:
```
.
├── R
│   ├── malignant
│   ├── microenvironment
│   ├── revisions
│   └── utils
├── README.md
├── data
├── genome-biology-methods.Rproj
├── output
├── poetry.lock
├── pyproject.toml
├── python
│   ├── chromatin_analysis.ipynb
│   ├── doubles_coefficient_analysis.ipynb
│   └── gl261_cinemaot.ipynb
├── renv.lock
├── shared_data
│   └── deseq_output
└── shell
```
- `R` contains R code
- `python` contains Python code
- `shell` contains shell scripts for CUT&TAG analysis
- `shared_data/deseq_output` contains comprehensive DESeq2 outputs for our perturb-seq data relative to non-targeting same treatment condition controls. This data was generated
using the R [`DElegate`](https://github.com/cancerbits/DElegate) package.

## Maintenance
Authored by Christopher Zou, Ashir Borah, and John Liu. Feel free to fork/ping us through issues for questions. 
