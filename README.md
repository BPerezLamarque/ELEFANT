# ELEFANT



This tutorial explains how to run ELEFANT.


<p align="center">
    <img src="https://github.com/BPerezLamarque/ELEFANT/blob/main/example/ELEFANT.png" width="600">
</p>

<p align="center">
    <b>Figure 1: ELEFANT: Evolution of LatEnt traits for Ancestral Network reconsTruction.</b> <small>First, given a bipartite interaction network (N), ELEFANT estimates the latent traits of extant species from both clades (L and R)  responsible of present-day interactions using the random dot graph product (RDPG) method. Each species is characterized by a vector of d latent trait values; here for clarity only the first latent trait value is represented with the color code. Second, for each guild, ELEFANT adds the extinct and unsampled lineages on the phylogenetic trees (colored in green) using data augmentation, which assumes constant rates of speciation (lambda), extinction (mu), and sampling at present (rho). Third, ELEFANT performs ancestral reconstructions of the latent traits on the phylogenetic trees of each guild using Brownian motions. Fourth, at any time in the past (e.g. 50 Mya), ELEFANT outputs ancestral interaction networks from the estimated ancestral latent trait values. Fifth, to assess whether the assumptions underlying ELEFANT are appropriate given a specific empirical network, a cross-validation procedure measures the ability of ELEFANT to reliably recover known interactions at present in this given system.</small> 
</p>


**Citation:** Benoît Perez-Lamarque, Jérémy Andréoletti, Baptiste Morillon, Orane Pion-Piola, Amaury Lambert, and Hélène Morlon, *Darwin’s Entangled Bank Through Deep Time*, bioRxiv, 2025, DOI: https://doi.org/10.1101/2025.10.08.681159



**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com




# Contents:
**[Installation](#installation)**\
**[Running ELEFANT](#running-elefant)**





# Installation:


The R script to run ELEFANT ("functions_ELEFANT.R") must be downloaded from the folder ['R'](https://github.com/BPerezLamarque/ELEFANT/tree/main/R/) and stored in your R working directory. 
In addition, the following packages must be installed:

```r
install.packages("phytools", "mvMORPH", "RPANDA", "bipartite", "igraph", "ggplot2", "ggpubr")

```



# Running ELEFANT:


First, you can set your working directory and **load the packages and functions**:

```r

setwd("YOUR_WORKING_DIRECTORY")


library(phytools)
library(RPANDA)
library(bipartite)
library(dplyr)
library(tidyr)
library(igraph)
library(reshape2)
library(ggplot2)
library(parallel)
library(ggpubr)

source("functions_ELEFANT.R")

```

Then, you can **load the example dataset** of the plant-Nymphalini interactions (from Braga et al., 2020) that can be downloaded from the folder ['example'](https://github.com/BPerezLamarque/ELEFANT/tree/main/example/): 

```r

# Open the phylogenetic trees:
tree_B <- read.tree("tree_Nymphalini.tre")
tree_A <- read.tree("tree_angiosperm_families.tre")
# Warning: Please note that the phylogenetic tree must be rooted, binary, and ultrametric to run ELEFANT. 

# Open the interaction network:
network <- read.table("network_Nymphalini_plants.csv", sep=";", header=TRUE)
# Warring: Each row corresponds to a taxon of clade B(Nymphalini here) and each column to a taxon of clade A (Angiosperm family)

```

Next, you need to **chose different parameters** before running ELEFANT:

```r

name <- "Nymphalini_plant" # select the name of the ELEFANT run

# Are interactions obligate? 
obligate_A=FALSE # Plants do not necessarily interact with herbivorous butterfly 
obligate_B=TRUE # Butterflies obligatorily interact with at least one plant species (they need to feed!)

# List of past ages for reconstructing the ancestral networks.
list_ages <- c(seq(0, 22, 2)) # Warning: it can not be older than the youngest MRCA of clades A and B. Here the Nymphalini clade is 22.3 Myr old. 

# Number of data augmentation for the network reconstructions
nb_recon <- 25 # must be at least 250; here 25 reconstructions is just to provide a fast-to-run example

# Indicate which clades to subsample for cross-validation (Step 5), and what proportion to use?
perc_cv_A=0.1
perc_cv_B=0.1

# Measure of global metrics
global_metrics=TRUE # whether you want to compute the global metrics of the ancestral networks (connectance, nestedness, modularity...). This step can be quite long! Set to FALSE if you are only interested in the ancestral interactions. 
null_model=TRUE # whether you want to compute the global metrics to null expectations (highly recommanded if you plan to interpret the significance and temporal trends of the global metrics).


```

*ELEFANT can be run* with the following function:


```r

results <- ELEFANT(name=name, nb_recon=nb_recon, list_ages=list_ages, 
                   tree_A=tree_A, tree_B=tree_B,
                   data_augmentation_A=TRUE, data_augmentation_B=TRUE, 
                   obligate_A=obligate_A,
                   obligate_B=obligate_B,
                   lambda_A=0.01, mu_A=0.001, rho_A=0.9, # Diversification rates used for the plants (values chosen only for the sake of example).
                   lambda_B=0.1, mu_B=0.01, rho_B=0.7, # Diversification rates used for the butterflies (values chosen only for the sake of example).
                   only_CV=FALSE, # whether only cross-validations (step 5) should be done (skipping steps 2-4)
                   perc_cv_A=perc_cv_A, perc_cv_B=perc_cv_B,
                   global_metrics=global_metrics,
                   null_model=null_model)

```

Finally, the followings *plots of ancestral networks and their associated global metrics* can be represented with the following formulas:

```r

plot_networks_ELEFANT(name, results)

plot_metrics_ELEFANT(name, results,
                     clade_A="clade_A",
                     clade_B="clade_B",
                     min_age=10)

```

All results are stored in a specific folder in your working directory. 
