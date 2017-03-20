Status: Travis CI [![Travis-CI Build Status](https://travis-ci.org/ipw012/RIVER.svg?branch=master)](https://travis-ci.org/ipw012/RIVER)

RIVER
=======

`RIVER` is an `R` package of a probabilistic modeling framework, __RIVER (RNA-Informed Variant Effect on Regulation)__ that jointly analyzes personal genome (WGS) and transcriptome data (RNA-seq) to estimate the probability that a variant has regulatory impact in that individual. It is based on a generative model assuming that genomic annotations, such as the location of a variant with respect to regulatory elements, determine the prior probability that variant is a _functional regulatory variant_, which is an unobserved variable. The functional regulatory variant status then influences whether nearby genes are likely to display outlier levels of gene expression in that person. __RIVER__ is trained in an unsupervised manner such that it does not require a labeled set of functional/non-functional variants; rather it derives its predictive power from identifying expression patterns that tend to coincide with particular rare variant annotations.

For more information about `RIVER`, check the vignettes.

For further details of a list of genomic annotations used for constructing features and how to generate the features and outlier status, please refer to our [submitted publication](http://biorxiv.org/content/early/2016/09/09/074443).

# Installation

Get most recent version of R (>= 3.3.2) from [CRAN](http://cran.r-project.org/).

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RIVER")
```
# Vignettes

The vignettes for this package can be viewed via [Bioconductor's website](http://www.bioconductor.org/packages/RIVER) (manual [backup](https://github.com/ipw012/RIVER)).

# Citation

Below is the citation output from using `citation('RIVER')` in R. __Please run this yourself to check for any updates on how to cite RIVER__.

To cite the `RIVER` package in publications use:

Xin Li, Yungil Kim, Emily K. Tsang, Joe R. Davis, Farhan N. Damani, Colby Chiang, Zachary Zappala, Benjamin J. Strober, Alexandra J. Scott, Andrea Ganna, Jason Merker, GTEx Consortium, Ira M. Hall, Alexis Battle, Stephen B. Montgomery (2016), "The impact of rare variation on gene expression across tissues", doi: https://doi.org/10.1101/074443, (URL: https://doi.org/10.1101/074443), <URL: http://biorxiv.org/content/early/2016/09/09/074443>.

A BibTeX entry for LaTeX users is

@article {Li074443,
	author = {Li, Xin and Kim, Yungil and Tsang, Emily K. and Davis, Joe R. and Damani, Farhan N. and Chiang, Colby and Zappala, Zachary and Strober, Benjamin J. and Scott, Alexandra J. and Ganna, Andrea and Merker, Jason and , and Hall, Ira M. and Battle, Alexis and Montgomery, Stephen B.},
	title = {The impact of rare variation on gene expression across tissues},
	year = {2016},
	doi = {10.1101/074443},
	publisher = {Cold Spring Harbor Labs Journals},
	journal = {bioRxiv},
	URL = {http://biorxiv.org/content/early/2016/09/09/074443}
}

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/) as well as Bioconductor's nightly build.
