# CryptosporidiumSlippage
A repository for the analyses of PCR slippage in _Cryptosporidium_ spp. gp60 amplicons

## Overview

A repository for `R` and `bash` scripts to accompany our manuscript on our observations on PCR slippage in _Cryptosporidium_ spp. gp60 amplicons using the `dada2` [tutorial](https://benjjneb.github.io/dada2/tutorial.html) as a basis, but with changes to key parameters.  Such changes are described below.

**Title**: "Quantifying replication slippage error in _Cryptosporidium_ metabarcoding studies"

**Authors**: Matt A Knox, Patrick J Biggs, Juan-Carlos Garcia-R and David T S Hayman.

This repository contains the following files and scripts:
* `01dataGeneration.R`
* `02bashCode.sh`
* `03dataPlotting.R`
* `Samples.txt`

## Requirements to run the scripts 

### environment
The scripts were run on a Ubuntu 20.04.6 LTS laptop with 32Gb RAM and an Intel® Core™ i7-7500U CPU @ 2.70GHz processor.

### R environment
The code was run in `R` (v. 4.3.1) from an RStudio environment (2023.06.0, Build 421) with the following required packages:
* dada2 (v. 1.26.0)
* phyloseq (v. 1.42.0)
* tidyverse (v. 1.3.2)

### Bash environment
In addition, the `bash` script between the two `R` scripts was run on in an Ubuntu 20.04 environment along with cd-hit for clustering and generating tabular output:
* [cdhit on the web](https://sites.google.com/view/cd-hit) and [cdhit in GitHub](https://github.com/weizhongli/cdhit) and associated scripts that is set to run from the path.


## Script description

### 01dataGeneration.R
The first part of the process uses `dada2` in R to generate the output.  All processes are initally run in a folder where the source data is, and the output can be written to there, or can be placed in a new folder elsewhere.  If that second option is taken, that folder can be made prior to data output, or inside `R`.  


### 02bashCode.sh
Once `01dataGeneration.R` is run, a piece of code to run to convert the sequences from dada2 into a text file for input back in the R environment.  This is run from the same folder where the dada2 results are.  It requires that cd-hit be accessible in the path, as  are all the extra scripts for working with cd-hit output.


### 03dataPlotting.R
This is code to generate Figure 1 from the manscript. 

### Samples.txt
This is a file used in plotting of the data in the generation of Figure 1. It has the following columns:

| Column | Purpose |
|--|--|
| Name | Sample name |
| WellNumber | The well numers the samples were run in |
| Description | A short description of the samples |
| Replicate | Technical replicate number |
| AltName | Alternative sample name |
| Group | PCR Group |
| Hominis | Proportion of _C. hominis_ in the sample mixture |
| Parvum | Proportion of _C. parvum_ in the sample mixture |
| HPMix | Ratio of _C. hominis_ to _C. parvum_ |

