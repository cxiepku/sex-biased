# sex-biased
evolution of sex-biased gene expression

This repository includes the R code and test data for the manuscript "Fast evolutionary turnover and overlapping variances of sex-biased gene expression patterns defy a simple binary classification of sexes".

The R script "analyze_sex.R" identifies and analyzes sex-biased genes, including the calculation of "IQR / Median" and "sex-bias index (SBI)". It takes a TPM matrix with gene annotation (example: "data/DOM_brain_TPM.tsv"), a sample annotation file (example: "data/DOM_brain_coldata.tsv") as input files. It also takes a few input parameters (see the script). It generates two sets of output files. One includes all the sex-biased and unbiased genes (example: "data/DOM_brain_sex.tsv"; in column "Type", "F" means female-biased, "M" means male-biased, and "U" means unbiased; column "Disp" is for "IQR / Median"), and the other includes SBI and scaled SBI values, and the density plot for the scaled SBI values for each sample (examples: "data/DOM_brain_SBI.tsv" and "data/DOM_brain_SBI.pdf"). There is also another group of examples for MUS brain in the folder "data".

The R script "analyze_turnover.R" identifies regulatory turnover genes between taxa. It takes a TPM matrix with gene annotation and a sample annotation file for each of the two taxa as input files. It also takes a few input parameters (see the script), and its algorithm is quite similar with the above sex-bias one. It generates one output file, including the information of expression (median TPM > 1 in at least one taxon, column "F_filtered" or "M_filtered") and turnover (column "F_turnover" or "M_turnover") for each of the female and male swaps (example: "data/DOM_MUS_brain_turnover.tsv").

The R scripts only require three R packages: matrixStats, dplyr, and ggplot2, and runs very fast. It has been tested under R 4.2.3, matrixStats 1.3.0, and dplyr 1.1.4, ggplot2 3.4.4, but the versions do not matter much.

The test data mentioned above are in the folder "data/", which are the mouse brain data of two taxa (see the manuscript for further information). The full data is in the supplemental materials of the manuscript, Edmond, or ENA (see the manuscript).
