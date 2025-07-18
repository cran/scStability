% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scStability.r
\name{scStability}
\alias{scStability}
\title{A user friendly wrapper function that runs the entire scRNA-seq stability workflow and shows statistics for each step}
\usage{
scStability(
  seurat_obj,
  n_runs = 100,
  dr_method = "umap",
  clust_method = "louvain",
  n_cores = 1,
  verbose = TRUE,
  print_plot = TRUE,
  seeds = NULL
)
}
\arguments{
\item{seurat_obj}{A Seurat object containing scRNA-seq data and a PCA}

\item{n_runs}{Number of DR embeddings and number of cluster assignments to be generated (< 250 recommended)}

\item{dr_method}{Method to use for dimension reduction, either "umap" or "tsne"}

\item{clust_method}{Algorithm used for clustering, either "louvain" or "leiden"}

\item{n_cores}{Number of CPU cores to use for parallelising functions}

\item{verbose}{Whether the function should print summary statistics as it calculates them}

\item{print_plot}{Whether the final medoid plot should be printed}

\item{seeds}{A set of seeds of length n_runs used for generating embeddings and clusters}
}
\value{
A list containing:
\item{mean_emb}{Data frame containing the mean embedding coordinates}
\item{mean_clust}{Vector of the mean cluster assignments}
\item{plot}{ggplot2 object with the medoid embedding plot and cluster assignments}
\item{embedding_stats}{List of embedding statistics}
\item{cluster_stats}{List of clustering statistics}
\item{seurat_object}{Seurat object now containing mean embeddings and mean clusters}
}
\description{
A wrapper function that runs all other stability analysis functions in order.
Statistics for each step are printed accordingly and a final DR and cluster plot is shown
which represents the medoid embeddings and cluster assignments that were generated.
}
