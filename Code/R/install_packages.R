# Set my personal .libPaths to install things correctly
# .libPaths(new="/usr/local/lib/R/site-library/3.5")


if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

# Data format package
BiocManager::install("SummarizedExperiment", version = "3.8")
BiocManager::install("SingleCellExperiment", version = "3.8")

# Basic Dimensionality reduction (PCA)
BiocManager::install("pcaMethods", version = "3.8")

# Single-cell RNA seq analysis package (Aaron Lun, Cambridge!)
BiocManager::install("scran", version = "3.8")

# Plot the results of PCA
#install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")

# Sparse PCA
install.packages("sparsepca")

# Advanced Dimensionality reduction (Isomap, LLE)
BiocManager::install("RDRToolbox", version = "3.8")

# Piping operation
install.packages("magrittr")
install.packages("tidyverse")

# Plotting
installed.packages("ggplot2")
install.packages("plotly")

# Pretty printing
install.packages("knitr")
install.packages("printr")

# Copy installed packages from personal folder to shared folder then delete from personal folder
# sudo cp -r ~/R/x86_64-pc-linux-gnu-library/3.5/* /usr/local/lib/R/site-library/3.5/
# rm -rf  ~/R/x86_64-pc-linux-gnu-library/3.5/*
# Check libpaths: .libPaths(), make sure it has "/usr/local/lib/R/site-library/3.5"



