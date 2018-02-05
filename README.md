Status: Travis CI [![Travis-CI Build Status](https://travis-ci.org/ipw012/RIVER.svg?branch=master)](https://travis-ci.org/ipw012/RIVER)

RIVER
=======

`RIVER` is an `R` package of a probabilistic modeling framework, __RIVER (RNA-Informed Variant Effect on Regulation)__ that jointly analyzes personal genome (WGS) and transcriptome data (RNA-seq) to estimate the probability that a variant has regulatory impact in that individual. It is based on a generative model assuming that genomic annotations, such as the location of a variant with respect to regulatory elements, determine the prior probability that variant is a _functional regulatory variant_, which is an unobserved variable. The functional regulatory variant status then influences whether nearby genes are likely to display outlier levels of gene expression in that person. __RIVER__ is trained in an unsupervised manner such that it does not require a labeled set of functional/non-functional variants; rather it derives its predictive power from identifying expression patterns that tend to coincide with particular rare variant annotations.

For more information about `RIVER`, check the vignettes.

For further details of a list of genomic annotations used for constructing features and how to generate the features and outlier status, please refer to our [submitted publication](http://biorxiv.org/content/early/2016/09/09/074443).

# Installation

Get most recent version of R (>= 3.4) from [CRAN](http://cran.r-project.org/).

```R
## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("RIVER")
```
# Vignettes

The vignettes for this package can be viewed via [Bioconductor's website](http://www.bioconductor.org/packages/RIVER) (manual [backup](https://github.com/ipw012/RIVER)).

# Citation

Below is the citation output from using `citation('RIVER')` in R. __Please run this yourself to check for any updates on how to cite RIVER__.

To cite the `RIVER` package in publications use:

Li, X , Kim, Y , Tsang, EK , Davis, JR , Damani, FN, Chiang, C, Hess, GT, Zappala, Z, Strober, BJ, Scott, AJ, Li, A, Ganna, A, Bassik, MC, Merker, J, GTEx Consortium, Hall, IM, Battle, A , Montgomery, SB , "The impact of rare variation on gene expression across tissues", Nature 550, 239-243 (2017), doi:10.1038/nature24267, <URL: https://www.nature.com/articles/nature24267>.

# Testing

Testing on Bioc-devel is feasible thanks to [R Travis](http://docs.travis-ci.com/user/languages/r/) as well as Bioconductor's nightly build.
