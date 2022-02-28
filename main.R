#Load Packages
library(DESeq2)
library(tidyverse)

#' Load a tsv located at specific location `filename` into a tibble
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (g x m) tibble with a 'gene' column followed by
#' sample names as column names. 
#' 
#' @note Column 'gene' should be first and the only column to contain strings.
#' Data in sample_name columns CANNOT be strings
#' 
#' @example `verse_counts <- read_data('verse_counts.tsv')`

read_data <- function(filename){
  verse_counts_tbl <- read_tsv(filename)
  head(verse_counts_tbl)
  return(verse_counts_tbl)
}

#' Filter out genes with zero variance
#'
#'
#' @param verse_counts tibble: a (g x m) tibble of raw read counts
#'
#' @return tibble: a (n x m) tibble of raw reads with genes that have 
#' zero variance across samples removed
#' 
#' @note (g >= n)
#' 
#' @example `filtered_counts <- filter_zero_var_genes(verse_counts)`

filter_zero_var_genes <- function(verse_counts) {
  verse_counts$var <- rowVars(as.matrix(verse_counts[,c(-1)]))
  pos <- which(verse_counts$var!=0)
  return(select(verse_counts, -var)[pos,])
}

#' Extract time point information from sample name
#'
#'
#' @param str string: sample name from count data.
#'
#' @return string: string character representing sample time point 
#'
#' @example `timepoint_from_sample("vAd_1")`
#' output:`"Ad"`

timepoint_from_sample <- function(str) {
  return(substr(str, 2, 3))
}


#' Grab sample replicate number from sample name
#'
#' @param str  string: sample name from count data.
#'
#' @return string: string character represent sample replicate number
#' 
#' @note you may choose to return numeric values instead of strings here
#'
#' @example `sample_replicate("vAd_1")`
#' output: `"1"`

sample_replicate <- function(str) {
  return(substr(str, 5, 5))
}

#' Generate sample-level metadata from sample names. 
#' 
#' Will include columns named "sample", "timepoint", and "replicate" that store 
#' sample names, sample time points, and sample replicate, respectively.
#'
#'
#' @param sample_names vector: character vector of length (_S_) consisting of sample
#' names from count data.
#'
#' @return tibble: a (_S_ x 3) tibble with column names "sample", 
#' "timepoint", and "replicate". "sample"holds sample_names; "timepoint" 
#' stores sample time points; and "replicate" stores sample replicate
#' 
#' @note _S_ < m
#'
#' @example `meta <- meta_info_from_labels(colnames(count_data)[colnames(count_data)!='gene'])`

meta_info_from_labels <- function(sample_names) {
  timepoint <- timepoint_from_sample(sample_names)
  replicate <- sample_replicate(sample_names)
  return(tibble(sample=sample_names, timepoint, replicate))
}

#' Calculate total read counts for each sample in a count data.
#'
#'
#' @param count_data tibble: a (n x m) tibble of raw read counts.
#'
#' @return named vector: numeric vector of read totals from each sample
#'
#' @examples `get_library_size(count_data)`

get_library_size <- function(count_data) {
  lib_size <- apply(count_data[-1], 2,sum)
  return(lib_size)
}

#' Normalize raw count data to counts per million WITH pseudocounts using the 
#' following formula:
#'     (count + 1) / ((sample_library_size/10^6) + 1)
#'
#'
#' @param count_data tibble: a (n x m) matrix of raw read counts.
#'
#' @return tibble: a (n x m) matrix with read count normalized to counts
#' per million
#'
#' @examples
#' `normalize_by_cpm(count_data)`

normalize_by_cpm <- function(count_data) {
  #3 a.
  cpm <- apply(count_data[,c(-1)],2, function(x) (x/sum(x))*1000000)
  cpm_tbl <- as.tibble(cpm)
  cpm_tbl$gene <- count_data$gene
  cpm_tbl <- relocate(cpm_tbl, gene, 1)
  return(cpm_tbl)
}

#' Normalize raw count data using DESeq2
#'
#'
#' @param count_data tibble: a (n x m) matrix of raw reads
#' @param meta_data tibble: sample-level information tibble containing time point,
#' sample, and replicate information.
#' @param design_formula formula: formula of comparision of interest
#'
#' @return tibble: a (n x m) tibble of DESeq2 normalized count data.
#'
#' @example ' `deseq_normalize(count_data, meta_data, ~ timepoint)`

deseq_normalize <- function(count_data, meta_data, design_formula) {
  # if I do design = design_formula, the test was giving me a warning, hence I did design = ~1
  # as suggested in the assignment sample report
  # dds <- DESeqDataSetFromMatrix(countData = count_data[-1], colData = meta_data, design = design_formula)
  dds <- DESeqDataSetFromMatrix(countData = count_data[-1], colData = meta_data, design = ~1)
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized=TRUE)
  norm_counts_tbl <- as.tibble(norm_counts)
  norm_counts_tbl$gene <- count_data$gene
  norm_counts_tbl <- relocate(norm_counts_tbl, gene, 1)
  return(norm_counts_tbl)
}

#' Perform and plot PCA using processed data.
#' 
#' PCA is performed over genes, and samples should be colored by time point.
#' Both `y` and `x` axis should have percent of explained variance included.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param meta tibble: sample-level meta information (_S_ x 3)
#' @param title string: title for plot
#'
#' @return ggplot: scatter plot showing each sample in the first two PCs. 
#'
#' @examples
#' `plot_pca(data, meta, "Raw Count PCA")`

plot_pca <- function(data, meta, title="") {
  t_data <- t(data)
  pca <- prcomp(t_data)
  pca_var <- pca$sdev**2 / sum(pca$sdev**2)
  x_label <- paste0("PC1: ",round(pca_var[1]*100, 2),"% variance")
  y_label <- paste0("PC2: ",round(pca_var[2]*100, 2),"% variance")
  pca_plot_data <- meta
  pca_plot_data$PC1 <- pca$x[,1]
  pca_plot_data$PC2 <- pca$x[,2]
  pca_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, col=timepoint)) +
    geom_point() +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(title)
  return(pca_plot)
}

#' Plot gene count distributions for each sample using boxplots.
#' 
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale the `y` axis to log10 values.
#' Default is FALSE, and y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: boxplot show gene count distributions for each sample
#'
#' @example `plot_sample_distributions(data, scale_y_axis=TRUE, title='Raw Count Distributions')`

plot_sample_distributions <- function(data, scale_y_axis=FALSE, title="") {
  long_data <- data %>% gather(samples, counts)
  if(scale_y_axis) {
    long_data <- filter(long_data, counts!=0)
  }
  gene_count_box_plot <- ggplot(long_data, aes(x=samples, y=counts, color=samples)) + 
    geom_boxplot() +
    ggtitle(title)
  if(scale_y_axis) {
    gene_count_box_plot <- gene_count_box_plot + scale_y_log10()
  }
  return(gene_count_box_plot)
}

#' Plot relationship between mean read counts and variability over all genes.
#'
#'
#' @param data tibble: a (n x _S_) data set
#' @param scale_y_axis boolean: whether to scale to y-axis to log10 values. Default
#' is false, and the y-axis will not be transformed.
#' @param title string: title to give the chart.
#'
#' @return ggplot: A scatter plot where the x-axis is the rank of gene ordered by mean
#' count over all samples, and the y-axis is the observed variance of the
#' given gene. Each dot should have their transparency increased. The scatter
#' plot should also be accompanied by a line representing the average mean and
#' variance values.
#'
#' @example `plot_variance_vs_mean(data, scale_y_axis=TRUE, title='variance vs mean (raw counts)')`

plot_variance_vs_mean <- function(data, scale_y_axis=FALSE, title="") {
  gene_means <- rowMeans(as.matrix(data))
  gene_vars <- rowVars(as.matrix(data))
  var_mean_data <- tibble(mean=gene_means, var=gene_vars)
  var_mean_data$rank <- rank(var_mean_data$mean)
  var_mean_plot <- ggplot(var_mean_data, aes(x=rank, y=var)) +
    geom_point(alpha=0.4) +
    geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
    xlab("Rank(Mean)") +
    ylab("Variance") +
    ggtitle(title)
  
  if (scale_y_axis) {
    var_mean_plot <- var_mean_plot + scale_y_log10()
  }
  return(var_mean_plot)
}
