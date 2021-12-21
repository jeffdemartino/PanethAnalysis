
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PanethAnalysis

**Author:** [Jeff DeMartino](https://github.com/jeffdemartino)
<a href="https://orcid.org/0000-0001-7366-4789" target="orcid.widget">
<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="16" height="16"/></a>

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
![GitHub tag (latest by
date)](https://img.shields.io/github/v/tag/jeffdemartino/PanethAnalysis)

<!-- badges: end -->

The `PanethAnalysis` package provides the data and wrapper functions
required to reproduce the analyses from our publication “Optimal human
intestinal organoid model reveals interleukin-22 induced Paneth cell
differentiation”.

## Installation

The following non-CRAN dependencies need to be installed prior to
PanethAnalysis:

``` r
# Bioconductor: DelayedMatrixStats, limma, destiny, slingshot
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DelayedMatrixStats", "limma", "destiny", "slingshot"))

# Github: SeuratWrappers
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github('satijalab/seurat-wrappers')
```

You can then install the development version of PanethAnalysis from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("jeffdemartino/PanethAnalysis", build_vignettes = TRUE)
```

Be sure to build the vignette upon installation (by invoking
`build_vignettes = TRUE`, as shown above) so the usage guide is made
available. Note that for the sake of convenience, the time/memory
intensive steps have been pre-computed (see
[here](https://ropensci.org/blog/2019/12/08/precompute-vignettes/) for
more details on how to replicate this).

As of this version (v0.1.0), all the required data sets are bundled in
this package (&gt;200MB). This may change in the future when the data is
made publicly available post-publication.

## Usage

The guide for using this package (i.e. reproducing the data analyses and
figure panels) is available as a package vignette, which can be accessed
by calling `browseVignettes("PanethAnalysis")` and selecting the HTML
file.

## Citation

XXX
