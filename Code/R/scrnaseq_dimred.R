# Perform various dimensionality reduction analyses on scRNA-seq data
setwd("~/Dropbox (Personal)/Documents/Applications/Cambridge2019_Gaurav/Cambridge_talk_20190306/Code/R")

# Load libraries
library(magrittr) # Pipe %>% operation for clean coding
library(SingleCellExperiment) # Data container (https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
library(pcaMethods) # Linear DimRed methods (PCA and extensions)
library(RDRToolbox) # Non-linear DimRed methods (LLE, Isomap)
library(ggbiplot)
library(plotly)

library(sparsepca)

#library(tidyverse)

# Read the data
data_zeisel <- readRDS("zeisel.rds")

print(data_zeisel)

# Gene statistics
rowData(data_zeisel)

# Cell statistics
colData(data_zeisel)

#For simplicity and speed, we work on a subset of 100 genes. To avoid ending up with only uninteresting genes, we extract the 100 genes with maximal variance in the log-transformed counts.
data_zeisel %>% 
  counts() %>%  # logcounts assay
  rowVars -> gene_count_variances # save row-wise variances

names(gene_count_variances) <- rownames(data_zeisel)

# Do the subsetting based on top 100 variances
data_zeisel_subset <- data_zeisel[
  names(
    (
      gene_count_variances %>% 
        sort(decreasing = TRUE)
      )[1:200]
    ),]
data_zeisel_subset




# Do basic PCA on raw subsetted count data
# results_pca <- pcaMethods::pca(
#   object = data_zeisel_subset %>% counts(),
#   method = "svd", 
#   nPcs = min(dim(data_zeisel_subset)), # Get all PCs
#   scale = "none",
#   center = TRUE, 
#   completeObs = FALSE, 
#   subset = NULL, 
#   cv = "none"
# )

results_princomp <- princomp(
  data_zeisel_subset %>% logcounts() %>% t()
)

results_prcomp <- prcomp(
  data_zeisel_subset %>% logcounts() %>% t()
)
results_prcomp_normalcounts <- prcomp(
  data_zeisel_subset %>% counts() %>% t()
)

summary(results_princomp)


results_spca <- sparsepca::spca(
  data_zeisel_subset %>% logcounts() %>% t(),
  alpha = 1e-2
)
rownames(results_spca$loadings) <- rownames(data_zeisel_subset)


results_spca_show <- results_princomp

for (list_elem in names(results_spca)){
  results_spca_show[[list_elem]] = results_spca[[list_elem]] 
}


# Plot PCA and Color by identified group (result of clustering)
P2 <- ggbiplot(results_prcomp_normalcounts, #results_spca_show,
               obs.scale = 1, 
               var.scale=1,
               ellipse=F,
               circle=F,
               varname.size=0.1,
               var.axes=F,
               groups=data_zeisel_subset$cell_type1, 
               alpha=0) 
P2$layers <- c(geom_point(aes(color=data_zeisel_subset$cell_type1), 
                          cex=1.0, alpha=0.3), P2$layers)

P2 %>% ggplotly()

# Get important genes
for (i in 1:3){
  cat("\n\nPC ", i, "\n\n")
  print(
    #results_prcomp$rotation[,i][order(abs(results_prcomp$rotation[,i]), decreasing = TRUE)][1:8]
    results_spca$loadings[,i][order(abs(results_spca$loadings[,i]), decreasing = TRUE)][1:20]
  )
}



# plot(results_pca)

set.seed(1234)
sample_size = 1500
sample_size = min(sample_size, dim(data_zeisel_subset)[2])
cell_sample <- sample(1:dim(data_zeisel_subset)[2], sample_size, replace=FALSE) # Sample a few cells randomly

# Run more advanced methods
if (sample_size <= 500){
results_isomap <- RDRToolbox::Isomap(
  data = (data_zeisel_subset[,cell_sample] %>% logcounts() %>% t()),
  dims = 1:10,
  k = 5,
  mod = FALSE,
  plotResiduals = TRUE,
  verbose = TRUE
  )
}

results_lle <- lapply(1:10, 
  function(x)RDRToolbox::LLE(
    data = (data_zeisel_subset[,cell_sample] %>% logcounts() %>% t()),
    dim = x,
    k = 8
)
)

results_prcomp_sample <- prcomp(
  (data_zeisel_subset[,cell_sample] %>% logcounts() %>% t())
)

# summary(results_prcomp_sample)
# results_prcomp_sample_residual_variance = (1-summary(results_prcomp_sample)$importance[3,]) #* sum(y=results_prcomp_sample$sdev^2)
# plot(x=1:length(results_prcomp_sample_residual_variance), y = results_prcomp_sample_residual_variance)

# "Residual variance" for PCA results (as defined weirdly in ISOMAP as 1 - cor(X_pw_dist, Y_pw_dist)^2)
X_pw_dist = as.matrix(dist(x=(data_zeisel_subset[,cell_sample] %>% logcounts() %>% t())))
resvar_pca = c()
for (i in 1:dim(results_prcomp_sample$x)[2]){
  y_PCA_pw_dist = as.matrix(dist(results_prcomp_sample$x[,1:i]))
  resvar_pca = c(resvar_pca, 
                 1 - cor(matrix(X_pw_dist, prod(dim(X_pw_dist)),1), matrix(y_PCA_pw_dist, prod(dim(X_pw_dist)), 1))^2
                 )  
}
plot(x=1:length(resvar_pca), y=resvar_pca)

# Residual variance for isomap
resvar_isomap = c()
for (i in 1:length(results_isomap)){
  y_isomap_pw_dist = as.matrix(dist(results_isomap[[i]]))
  resvar_isomap = c(resvar_isomap, 
                 1 - cor(matrix(X_pw_dist, prod(dim(X_pw_dist)),1), matrix(y_isomap_pw_dist, prod(dim(X_pw_dist)), 1))^2
  )  
}

# Residual variance for lle
resvar_lle = c()
for (i in 1:length(results_isomap)){
  y_lle_pw_dist = as.matrix(dist(results_lle[[i]]))
  resvar_lle = c(resvar_lle, 
                    1 - cor(matrix(X_pw_dist, prod(dim(X_pw_dist)),1), matrix(y_lle_pw_dist, prod(dim(X_pw_dist)), 1))^2
  )  
}

results_isomap_show <-  prcomp(
  (data_zeisel_subset[,cell_sample] %>% logcounts() %>% t()),
  rank.=2
)
results_isomap_show$x <- results_isomap$dim2
results_isomap_show$x <- results_lle[[2]]

# Plot PCA and Color by identified group (result of clustering)
P2 <- ggbiplot(results_isomap_show, #results_prcomp_sample,
               obs.scale = 1, 
               var.scale=1,
               ellipse=F,
               circle=F,
               varname.size=0.1,
               var.axes=F,
               groups=(data_zeisel_subset[,cell_sample])$cell_type1, 
               alpha=0) 
P2$layers <- c(geom_point(aes(color=(data_zeisel_subset[,cell_sample])$cell_type1), cex=0.5), P2$layers)

P2 %>% ggplotly()


#results_isomap$dim1

RDRToolbox::plotDR(results_isomap$dim2)
