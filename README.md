
<!-- README.md is generated from README.Rmd. Please edit that file -->

# se.norway.shoredate.14c

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh///master?urlpath=rstudio)

This repository contains the data and code for our paper:

> Authors, (YYYY). *Comparing summed probability distributions of
> shoreline and radiocarbon dates from the Mesolithic Skagerrak coast of
> Norway*. Name of journal/book <https://doi.org/xxx/xxx>

Our pre-print is online here:

> Roalkvam, Isak <a href="https://orcid.org/0000-0001-6974-1374">
> <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
> and Solheim, Steinar <a href="https://orcid.org/0000-0001-8293-8147">
> <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
> (2023). *Comparing summed probability distributions of shoreline and
> radiocarbon dates from the Mesolithic Skagerrak coast of Norway*. Name
> of journal/book, Accessed 24 Sep 2023. Online at
> <https://doi.org/xxx/xxx>

### How to cite

Please cite this compendium as:

> Roalkvam, Isak <a href="https://orcid.org/0000-0001-6974-1374">
> <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
> and Solheim, Steinar <a href="https://orcid.org/0000-0001-8293-8147">
> <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
> (2023). *Compendium of R code and data for Comparing summed
> probability distributions of shoreline and radiocarbon dates from the
> Mesolithic Skagerrak coast of Norway*. Accessed 24 Sep 2023. Online at
> <https://doi.org/10.5281/zenodo.8373857>

## Contents

The **R** directory contains scripts for performing the analyses and
create the figures used in the paper.

The **analysis** directory contains:

- [:file_folder: paper](/analysis/paper): R Markdown source document for
  manuscript. It also has a rendered version of the manuscript,
  `paper.pdf`, suitable for reading.
- [:file_folder: data](/analysis/data): Data used in the analysis.
- [:file_folder: figures](/analysis/figures): Plots and other
  illustrations

## How to run in your browser or download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the compendium as a zip from from this URL:
[master.zip](/archive/master.zip). After unzipping: - open the `.Rproj`
file in RStudio - run `devtools::install()` to ensure you have the
packages this analysis depends on (also listed in the
[DESCRIPTION](/DESCRIPTION) file). - finally, open
`analysis/paper/paper.Rmd` and knit to produce the `paper.docx`, or run
`rmarkdown::render("analysis/paper/paper.Rmd")` in the R console

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.
