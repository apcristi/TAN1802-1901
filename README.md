# TAN1802-1901
Files:

"Stations surface" Excel file includes metadata for the surface depth (10 to 25m) of all stations.

Scripts:

1.- Dada2
Pipeline adapted from https://vaulot.github.io/tutorials/R_dada2_tutorial.html

Dada2 prerequisities:
- Install R libraries library(dada2) # Must use version >= 1.12 library(Biostrings) library(ShortRead) library(stringr) library(ggplot2) library(dplyr) library(tidyr) library(tibble) library(readr) library(purrr) library("optparse")
- The dada2 algorithm requires primers to be removed prior to processing. Using dada2 there are 2 possibilities:
  1) Remove by sequence, but dada2 does not allow for ambiguities
  2) Remove by position, which is not a problem for Illumina sequences but is a problem for 454
For complex situation we recommend to use cutadapt to remove the primers : (https://github.com/marcelm/cutadapt/releases)
Note: to install cutadapt on macOS, cannot use github link above. Need to install with Conda (https://cutadapt.readthedocs.io/en/stable/installation.html). For mac I ran all the code to install cutadapt on terminal.
  1) Download latest miniconda here (https://docs.conda.io/en/latest/miniconda.html#installing).
  2) Follow this guide to install conda on macOS (https://conda.io/projects/conda/en/latest/user-guide/install/macos.html).
  3) Add bioconda channel (https://bioconda.github.io/user/install.html).
  4) install cutadapt.
     
2.- Phyloseq
- Pipeline adapted from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
- Input files (excel files)
  1) Metadata
  2) Otu table_1901-1802
  3) Tax table_1901-1802
     
3.- MicroViz 
- Pipeline adapted from https://david-barnett.github.io/microViz/

4.- MixOmics
- Pipeline adapted from http://mixomics.org/case-studies/splsda-srbct-case-study/
- Input files (excel files)
  1) 1802_Ps_correlation_OTU (ASV table TAN1802)
  2) TAN1802 (Metadata TAN1802)
  3) 1901_Ps_correlation_OTU (ASV table TAN1901)
  4) TAN1901 (Metadata TAN1901)
