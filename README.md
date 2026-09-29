# ChrAccR <img src="man/figures/chraccr_logo.png" align="right" height="96"/>


Welcome to `ChrAccR`, an R package that provides tools for the comprehensive analysis chromatin accessibility data. The package implements methods for data quality control, exploratory analyses (including unsupervised methods for dimension reduction, clustering and quantifying transcription factor activities) and the identification and characterization of differentially accessible regions. It can be used for the analysis of large bulk datasets comprising hundreds of samples as well as for single-cell datasets with 10s to 100s of thousands of cells. 

Requiring only a limited set of R commands, ChrAccR generates analysis reports that can be interactively explored, facilitate a comprehensive analysis overview of a given dataset and are easily shared with collaborators. The package is therefore particularly useful for users with limited bioinformatic expertise, researchers new to chromatin analysis or institutional core facilities providing ATAC-seq as a service. Additionally, the package provides numerous utility functions for custom R scripting that allow more in-depth analyses of chromatin accessibility datasets.

## Installation

To install `ChrAccR` and its dependencies, use the `devtools` installation routine:

```r
# install devtools if not previously installed
if (!is.element('devtools', installed.packages()[,"Package"])) install.packages('devtools')

# install dependencies
devtools::install_github("demuellae/muLogR")
devtools::install_github("demuellae/muRtools")
devtools::install_github("demuellae/muReportR")

# install ChrAccR
devtools::install_github("EpigenomeInformatics/ChrAccR", dependencies=TRUE)
```

## Getting started

The `ChrAccR` [vignette](https://epigenomeinformatics.github.io/ChrAccR/articles/overview.html) provides a most excellent starting point to get familiar with the package.

## How to cite

If you use `ChrAccR` in your work, please cite it:

> Mueller F, Gunduz IB (2026). ChrAccR: Analyzing chromatin accessibility data in R. R package version 0.9.26. https://github.com/EpigenomeInformatics/ChrAccR

```bibtex
@Manual{ChrAccR,
  title = {ChrAccR: Analyzing chromatin accessibility data in R},
  author = {Fabian Mueller and Irem B. Gunduz},
  year = {2026},
  note = {R package version 0.9.26},
  url = {https://github.com/EpigenomeInformatics/ChrAccR},
}
```

The entry for your installed version is available from R:

```r
print(citation("ChrAccR"), bibtex = TRUE)
```
