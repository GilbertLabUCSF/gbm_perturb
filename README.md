# GBM-Perturb

This is the analysis code for the gbm-perturb project. We're doing preliminary
set up at the moment, and some practices may change over time.

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

## Directory Structure