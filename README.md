# Search of active modules in protein-protein interaction networks

## Tasks of the project

* Study working of existing algorithms
* Learn to generate differential expression data to check active module searching algorithm developed in laboratory of bioinformatics (IFMO)
* Use active module searching algorithm on real data and evaluate its working in biological terms

## System requirements

Local computer with RStudio and following libraries:
BioNet, DLBCL, ggplot2, Hmisc, sqldf, igraph, dplyr, tidyr, MASS, affy, limma, oligo, pd.mogene.2.0.st,
genefilter, mouse4302.db, annotate, annaffy, knitr, mwcsr, mcmcRanking, rlist. 
Then, you must download necessary functions using this string:

```source("https://raw.githubusercontent.com/assaron/r-utils/master/R/exprs.R")```

## Result examples

![Active module](https://github.com/flowdel/diff_gene_expression/blob/master/ESC_OKSM_mcmc_last.pdf)

As a result you can obtain samples with average comparison p-values distributed as Beta-uniform, and this BU distribution shaped by 'a' value, which is defined by yourself.

**Simulated sample with Beta-uniform distributed p-values with 'a' value = 0.5**
![Beta-uniform distributed p-values with 'a' value = 0.5 ](https://github.com/flowdel/diff_gene_expression/blob/master/0_5_a.png)
**Simulated sample with Beta-uniform distributed p-values with 'a' value = 0.8**
![Beta-uniform distributed p-values with 'a' value = 0.8](https://github.com/flowdel/diff_gene_expression/blob/master/0_8_a.png)

## Databases

1. mouse4302 database for expression set annotation;
2. [InWeb-IM](https://omictools.com/inweb-inbiomap-tool) database of protein-protein ineractions.
