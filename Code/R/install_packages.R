if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Data format package
BiocManager::install("SingleCellExperiment", version = "3.8")

# Basic Dimensionality reduction (PCA)
BiocManager::install("pcaMethods", version = "3.8")

# Plot the results of PCA
library(devtools)
install_github("vqv/ggbiplot")

# Sparse PCA
install.packages("sparsepca")

# Advanced Dimensionality reduction (Isomap, LLE)
BiocManager::install("RDRToolbox", version = "3.8")

# Piping operation
install.packages("magrittr")




