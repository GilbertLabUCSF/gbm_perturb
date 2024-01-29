# GBM-Perturb

This repository contains anlaysis in Python and R for the output of _in vivo_ perturb-seq
experiments. Here's he directory structure:
```
.
├── cinema_ot/
│   └── requirements.txt
└── deseq_gsea_lda/
    ├── R/
    │   ├── malignant/
    │   └── microenvironment/
    ├── example_data/
    │   ├── malignant/
    │   └── microenvironment/
    ├── .Rprofile
    └── renv.lock
```
For dependency independence, the Python and R components are separated. Specifically, `cinema_ot`, for example,
contains Python code for CINEMA-OT analysis. The remaining analysis is in the `deseq_gsea_lda` directory
as is written in R.

## Setup

### `deseq_gsea_lda`

We use [`renv`](https://rstudio.github.io/renv/articles/renv.html) 
to manage dependencies and install libraries. Note, however, that `renv` is not 
always able to manage packages that have not been released to CRAN, 
Bioconductor, or other distribution platforms. `renv` also cannot manage non-R
packages, as it is run from within R. If you're working from RStudio, you may
need to modify your package discovery paths to include the path to `renv`'s
library of packages.

To get set up, first ensure that your working directory is the `deseq_gsea_lda` directory.

Then, install renv if you don't have it installed already:
```install.packages("renv")```

Next, initialize renv:
```renv::init()```

Restore the project dependencies from the `renv.lock` file:
```renv::restore()```

### `cinema_ot`

`requirements.txt` contains a dump of all the dependencies you might possibly need to run
the Jupyter notebook in `cinema_ot`. If you have pip, the easiest way to install them is 
to run the following command from the `cinema_ot` directory:
```
pip install -r requirements.txt
```

### Installing new packages

The `renv.lock` file is a snapshot of our R dependency world. When you add
additional dependencies, be sure to update it using the following command:
```renv::snapshot()```

If at any time you would like to return to the current state of `renv.lock`,
you can call the `renv::restore()` method from above.

## Maintenance

Code here was authored by Christopher Zou, Ashir Borah, and John
Liu. It is in a pretty scientifically useful, but pretty low standard
from a software engineering perspective state. Please feel free to fork or 
make PRs as you see fit in order to make things easier to work with!
