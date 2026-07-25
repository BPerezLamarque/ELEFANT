# ELEFANT



This tutorial explains how to run ELEFANT for reconstructing ancestral interaction network from present-day data.


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
**[Analyzing outputs of ELEFANT](#analyzing-outputs-of-elefant)**




# Installation:


The R script to run ELEFANT ("functions_ELEFANT.R") must be downloaded from the folder ['R'](https://github.com/BPerezLamarque/ELEFANT/tree/main/R/) and stored in your R working directory. 
In addition, the following packages must be installed:

```r
install.packages("phytools", "mvMORPH", "RPANDA", "bipartite", "igraph", "ggplot2", "purrr", "ggpubr")

```



# Running ELEFANT:


First, you can set your working directory and **load the required packages and functions**.
ELEFANT has been successfully tested on Linux and MacOS using R versions >4.2. 
The running time of ELEFANT will depend on the size of the datasets (it can take from a few minutes to several hours). To save time, the computation of the global metrics of the ancestral networks (the longest part) can be skipped. 

```r

setwd("YOUR_WORKING_DIRECTORY")


library(phytools)
library(RPANDA)
library(bipartite)
library(purrr)
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
# Load the phylogenetic trees:
tree_B <- read.tree("tree_Nymphalini.tre")
tree_A <- read.tree("tree_angiosperm_families.tre")

# Warning: the phylogenetic trees must be rooted, binary, and ultrametric to run ELEFANT.

# Load the interaction network:
network <- read.table("network_Nymphalini_plants.csv", sep = ";", header = TRUE)

# Warning: the interaction network must be provided as a binary interaction matrix (0 = no interaction, 1 = interaction). 
# Each row corresponds to a taxon from clade B (here, Nymphalini butterflies), and each column corresponds to a taxon from clade A (here, angiosperm families).

```

Next, you need to **choose the different parameters** before running ELEFANT:

```r
# Name of the ELEFANT run
name <- "Nymphalini_plant"

# Are interactions obligate?
obligate_A <- FALSE  # Plant families do not necessarily interact with Nymphalini butterflies.
obligate_B <- TRUE   # Butterflies obligatorily interact with at least one host plant species.

# Ages (in Myr) at which ancestral networks will be reconstructed
list_ages <- seq(0, 22, 2)

# Warning: the oldest age cannot exceed the age of the youngest MRCA of the two clades. Here, the Nymphalini clade is 22.3 Myr old.

# Number of ancestral network reconstructions
nb_recon <- 25

# Warning: at least 250 reconstructions are recommended for empirical analyses. Here, we use only 25 reconstructions to provide a fast-running example.

# Proportion of species to subsample during cross-validation (step 5)
perc_cv_A <- 0.1
perc_cv_B <- 0.1

# Compute global network metrics?
global_metrics <- TRUE

# If TRUE, connectance, nestedness, modularity, and other network metrics are computed for each ancestral network. This step can be computationally intensive. Set to FALSE if you are only interested in reconstructing ancestral interactions.

# Compare network metrics with null expectations?
null_model <- TRUE

# Highly recommended if you plan to interpret the significance or temporal dynamics of the reconstructed network metrics.
```

<br> <br>

**ELEFANT can then be run** using the following function:

```r
results <- ELEFANT(
  name = name,
  network = network,
  tree_A = tree_A, tree_B = tree_B,
  nb_recon = nb_recon,
  list_ages = list_ages,
  threshold = "Youden",
  only_CV = FALSE,   # If TRUE, only performs cross-validation (step 5), skipping steps 2-4.
  perc_cv_A = perc_cv_A, perc_cv_B = perc_cv_B,
  obligate_A = obligate_A, obligate_B = obligate_B,
  data_augmentation_A = TRUE, data_augmentation_B = TRUE,
  lambda_A = 0.01, mu_A = 0.001, rho_A = 0.9,  # Example diversification parameters for the plants.
  lambda_B = 0.10, mu_B = 0.010, rho_B = 0.7,  # Example diversification parameters for the butterflies.
  global_metrics = global_metrics,
  null_model = null_model
)
```

<br> <br>

Here are more details about the *arguments of the `ELEFANT()` function*:

| Argument         | Description                                                                                                                                                                                                                                                     |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`           | Name of the ELEFANT run.                                                                                                                                                                                                                                        |
| `network`        | A bipartite interaction matrix with species from guild A in columns and species from guild B in rows. The network must be binary (0 = no interaction, 1 = interaction).                                                                                                                                                           |
| `tree_A`         | Phylogenetic tree of guild A (corresponding to the columns of `network`). The phylogenetic tree must be rooted, binary, and ultrametric.                                                                                                                                                                                     |
| `tree_B`         | Phylogenetic tree of guild B (corresponding to the rows of `network`). The phylogenetic tree must be rooted, binary, and ultrametric.                                                                                                                                                                                     |
| `nb_recon`       | Number of ancestral network reconstructions to perform. We recommend using at least 250 reconstructions.                                                                                                                                                        |
| `list_ages`      | Vector of ages (in the same time units as the phylogenies) at which ancestral networks are reconstructed. The sequence may start at 0 (the present) and cannot extend beyond the age of the youngest MRCA of the two clades.                                    |
| `only_CV`        | If `TRUE`, only performs the cross-validation step (step 5), skipping the more computationally intensive reconstruction steps (steps 2-4).                                                                                                                      |
| `perc_cv_A`      | Proportion of species from clade A to subsample during cross-validation (step 5). A value of `0` indicates no subsampling for clade A, whereas `0.1` indicates that 10% of the species are subsampled.                                                                      |
| `perc_cv_B`      | Proportion of species from clade B to subsample during cross-validation (step 5). A value of `0` indicates no subsampling for clade B, whereas `0.1` indicates that 10% of the species are subsampled.                                                                      |
| `obligate_A`     | Whether species in clade A are obligate interactors (i.e. each species is constrained to interact with at least one partner).                                                                                                                                   |
| `obligate_B`     | Whether species in clade B are obligate interactors (i.e. each species is constrained to interact with at least one partner).                                                                                                                                   |
| `global_metrics` | If `TRUE`, computes global metrics of the reconstructed ancestral networks (e.g. connectance, nestedness, modularity). This step can be computationally intensive, so set it to `FALSE` if you are only interested in reconstructing ancestral interactions (and not the ancestral network structures). |
| `null_model`     | If `TRUE`, compares the global network metrics with null expectations. This option is highly recommended when interpreting the significance or temporal dynamics of network structure.                                                                          |
| `seed`           | Random seed used to ensure reproducibility of the analyses.                                                                                                                                                                                                     |
| `nb_cores`       | Number of CPU cores used to run the analyses (default = 1).                                                                                                                                                                                                     |
| `path_results`   | Directory where the results will be saved (default = the current working directory).                                                                                                                                                                            |

                
<br><br>

Then, the following arguments are *specific to the data augmentation procedure* for clade A (the corresponding arguments for clade B are analogous):

| Argument              | Description                                                                                                                                                                                                                                                                                                                                                                  |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `data_augmentation_A` | Specifies whether data augmentation should be performed for clade A. If `TRUE`, two options are available: (i) data augmentation using a constant-rate birth–death model, in which case `lambda_A`, `mu_A`, and `rho_A` must be provided; or (ii) supplying a list of pre-generated list of augmented trees obtained using a more complex diversification model (e.g. ClaDS or BDD). |
| `lambda_A`            | Speciation rate of clade A.                                                                                                                                                                                                                                                                                                                                                  |
| `mu_A`                | Extinction rate of clade A.                                                                                                                                                                                                                                                                                                                                                  |
| `rho_A`               | Sampling fraction of clade A, calculated as the number of species present in `tree_A` divided by the estimated total number of extant species in the clade.                                                                                                                                                                                                                  |
| `treesDA_A`           | A `multiPhylo` object containing pre-generated augmented phylogenetic trees for clade A, obtained using an alternative diversification model.                                                                                                                                                                                                                                |

<br> <br>


# Analyzing outputs of ELEFANT:

Finally, the followings **plots of ancestral networks and their associated global metrics** can be represented with the following formulas:

```r

plot_networks_ELEFANT(name, results)

plot_metrics_ELEFANT(name, results,
                     clade_A="Plants",
                     clade_B="Nymphalini",
                     min_age=10)

```

All results are stored in a specific folder in your working directory. 



<p align="center">
    <img src="https://github.com/BPerezLamarque/ELEFANT/blob/main/example/output_Nymphalini.png" width="600">
</p>

<p align="center">
    <b>Figure 2: Outputs of ELEFANT:</b> A. Changes of global metrics through time in the ancestral network. B. Consensus ancestral network of the plant-Nymphalini network 10 Myr ago.
</p>


