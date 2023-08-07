# GBM-Perturb

This is the analysis code for the gbm-perturb project.

## Setup

gbm-perturb uses [`renv`](https://rstudio.github.io/renv/articles/renv.html) 
to manage dependencies and install libraries. Note, however, that `renv` is not 
always able to manage packages that have not been released to CRAN, 
Bioconductor, or other distribution platforms. `renv` also cannot manage non-R
packages, as it is run from within R. If you're working from RStudio, you may
need to modify your package discovery paths to include the path to `renv`'s
library of packages.

To get set up, first ensure that your working directory is the project directory.

Then, install renv if you don't have it installed already:
```install.packages("renv")```

Next, initialize renv:
```renv::init()```

Restore the project dependencies from the `renv.lock` file:
```renv::restore()```

### Installing new packages

The `renv.lock` file is a snapshot of our R dependency world. When you add
additional dependencies, be sure to update it using the following command:
```renv::snapshot()```

If at any time you would like to return to the current state of `renv.lock`,
you can call the `renv::restore()` method from above.

## Directory Structure and File Guide

Current as of 8/7/2023
```
└── R/
    ├── utils.R
    ├── malignant/
    │   ├── nmf_full_noRT.R
    │   ├── nmf_full_RT.R
    │   ├── degenes.R
    │   ├── deseq.R
    │   ├── preprocess_malignant.R
    │   └── nmf_scripts/
    │       ├── analyze_nmf_RT.R
    │       └── run_nmf_RT.R
    └── microenvironment/
        ├── sb28_analysis.R
        └── preprocess_microenvironment.R
```
- `utils.R` contains useful functions for manipulating outputs
- `malignant/` houses all analysis for the GL261 data.
  - `preprocess_malignant.R`, `deseq.R`, `degenes.R`, and the two NMF files
    provide the bulk (in order) of our analysis.
  - The two files in `nmf_scripts/` are provided for ease of use when running
    NMF. You may find it easier to split up the running and visualizing steps.
- `microenvironment/` houses all analysis for the SB28 data.
  - `preprocess_microenvironment.R` includes data for setting up the Seurat
    object.
  - `sb28_analysis.R` includes deseq, DE gene extraction, and NMF

## Maintenance

The gbm-perturb project was written by Christopher Zou, Ashir Borah, and John
Liu. It is in a pretty scientifically useful, but pretty low standard
from a software engineering perspective state. Please feel free to fork or 
make PRs as you see fit in order to make things easier to work with!