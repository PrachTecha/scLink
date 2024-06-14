# scLink

## single-cell Link community detection (scLink)

This package is a wrapper around `linkcomm` package in R to perform the link community detection on single-cell omics data in Python. Additionally, this package is intended to be use as an extension of `Scanpy` package.

## Installation
It is a good practice to install this package in a clean environment with `conda`:

```sh 
    conda create -n scLink python=3.11 scanpy r-base r-essentials rpy2 -y
    conda activate scLink
    R -e 'install.packages("linkcomm",repos = "http://cran.us.r-project.org")'
```

Then you can install the package by:

```sh  
    pip install git+https://github.com/prachsk/scLink
```

to avoid any possible crashes due to rpy2 not finding the R install on conda, run the following import command:

```python
    import os, sys
    os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
    import scLink
```

## Example
For example usage of scLink please see `test.ipynb`.

## Citation
Kalinka, A.T. and Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. Bioinformatics 27 (14), 2011-2012. doi:10.1093/bioinformatics/btr311

    @Article{,
        title = {linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type},
        author = {Alex T Kalinka and Pavel Tomancak},
        year = {2011},
        journal = {Bioinformatics},
        volume = {27},
        pages = {2011-2012},
        doi = {10.1093/bioinformatics/btr311},
        }