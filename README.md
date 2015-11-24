This repository accompanies the paper:
Accumulation tests for FDR control in ordered hypothesis testing.
Ang Li and Rina Foygel Barber. arXiv:1505.07352

This paper proposes a family of accumulation tests for ordered hypothesis testing, where we seek to test a ranked list of n hypotheses, H1, ..., Hn, accompanied by p-values p1, ..., pn. The main idea is to estimate the false discovery proportion among the top k hypotheses in the ranked list, for each k=1,...,n, to determine an adaptive cutoff for the list. We aim to find a cutoff k that enables us to discover as many true signals as possible while controlling the false discovery rate. For more information, see the paper. 

Potential applications of this project are: 1. Allow us to leverage prior experiment and testing results to boost discoveries in new experimental data; 2. Provide a high dimensional linear model selection approach that controls the false detection of irrelevant predictors.

Here we give code in the R language providing an implementation of the accumulation test method, and reproducing the empirical results in the paper.
Description of files in this repository:

accumulation_test_functions.R -- R script to implement the accumulation test methods

ShortDemo.R, ShortDemo.pdf -- Here we create a sequence of hypotheses and p-values, and running accumulation test methods to test these hypotheses.  These two files are the R scripts and results (visualizations).

simulation.R, simulation.pdf -- Here we reproduce the simulated data results in the paper. These two files are the R scripts and results (visualizations).

gene_dosage_experiment.R, gene_dosage_experiment.pdf -- Here we reproduce the gene dosage data results in the paper. These two files are the R scripts and results (visualizations).
