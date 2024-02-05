
# ZINQ_L

<!-- badges: start -->
<!-- badges: end -->

This repo provides the code and example data to run ZINQ-L. 

## Overview
We propose a zero-inflated two-part quantile regression approach for differential abundance (DA) analysis of longitudinal microbiome data (ZINQ-L)

## Instructions
To run the code, you need to clone the repo.
```
git clone https://github.com/stephlee3/ZINQ_L.git
cd ZINQ_L
```

Install and load all the necessary packages.

```{r}
library(tidyverse)
library(phyloseq)
library(GUniFrac)
library(quantreg)
library(parallel)
library(lmerTest)
library(lme4)
library(CompQuadForm)
library(expm)
library(FactoMineR)
library(Matrix)
library(MASS)
library(LDM)
library(MicrobiomeStat)
library(Maaslin2)
library(NBZIMM)
```

Now you can run the example code in `example.R`. A matrix of P-values (each row: a replicate of representative taxa, each column: a statistical test methodology) will be generated on the simulated data. 

## Basic Structure
This repo is organized as follows

* `Code/`
  * `util.R`: contains all the statistical tests.
  * `README.md`: description of tests. 
  
* `example.R`: the example code to run the tests on simulated data. To run the sim, please contact the author to access the required data. 


## Contact Information
For any additional questions, please contact Runzhe Li(rli51@jhmi.edu).




