# cutandrun_SRSF2_Wagner2023

This repository contains a script to process cutandrun data as was done in Wagner et al. 2023. It takes as input BAM files from a cutandrun experiment and computes the pausing index for a list of genes, defined as the number of reads aligned to the TSS (+/-250bp) divided by the number of reads aligned to the gene body. In addition, it computes pvalues for the difference in pausing indices between two conditions.

Dependencies:

* numpy
* pandas
* scipy
* statsmodels
* matplotlib
* deeptools
* pysam
