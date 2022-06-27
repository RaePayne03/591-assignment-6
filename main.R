library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('reshape2')
# library('fgsea')


#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param subset(list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data
#' @export
#' 
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))

make_se <- function(counts_csv, metafile_csv, subset) {
  
  counts <- read.delim('verse_counts.tsv', sep ='\t', row.names ='gene') %>%
    data.matrix()
  
  meta_data <- readr::read_csv(metafile_csv) %>%
    dplyr::select(samplename, timepoint) %>%
    #filter for only vp0 and vAd
    dplyr::filter(timepoint %in% subset) %>%  
    #access components of factor using levels
    dplyr::mutate(timepoint = factor(timepoint, levels=subset)) 
  
  #pull returns a single column as a vector, select returns one or more columns as a data.frame -- doesnt work here.
  sample_names <- pull(meta_data, samplename) 
  new_subset <- counts [, sample_names]
  
  sum_exp <- SummarizedExperiment(assays = list(counts = new_subset), colData = meta_data)
  metadata(sum_exp)$model <- counts ~ timepoint
  
  return(sum_exp)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
#' 

return_deseq_res <- function(sum_exp, design) {
  
  #dds object returned by DESeq2
  dds <- DESeq2::DESeqDataSet(sum_exp,design)
  dds <- DESeq2::DESeq(dds)
  #convert result into data frame 
  out <- results(dds) %>% as.data.frame()
  #format output into list
  return(list('result' = out, 'dds' = dds))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq_result, padj_threshold){
  labeled_deseq <- deseq_result %>%
    as_tibble(rownames='genes') %>%
    #label logfold2change based on padj value
    dplyr::mutate(volc_plot_status = case_when(log2FoldChange < 0 & padj < padj_threshold ~ 'DOWN', 
                                               log2FoldChange > 0 & padj < padj_threshold ~ 'UP',
                                               TRUE ~ 'NS'))
  #reorder columns in tibble to match those in report
  dplyr::select(labeled_deseq,genes, volc_plot_status, log2FoldChange, padj, baseMean, lfcSE, stat,pvalue) %>%
    arrange(padj) %>%
    return()
} 


#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_plot <- labeled_results %>% 
    ggplot2::ggplot(aes(pvalue)) + 
    geom_histogram(bins=50, color='black', fill='lightblue') + 
    theme_light() + 
    ggtitle('Distribution of unadjusted p-values obtained from DE analysis (vP0 vs. vAd)') %>%
    return()
}


#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)

plot_log2fc <- function(labeled_results, padj_threshold) {
  
  logfc_plot <- labeled_results %>% 
    dplyr::filter(padj < padj_threshold) %>% 
    ggplot2::ggplot(aes(log2FoldChange)) + 
    geom_histogram(bins=100, color='black', fill='light blue') + 
    theme_light() + 
    ggtitle('Distribution of Log2FoldChanges for DE Genes (vP0 vs. vAd)') %>%
    return()
}



#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  #extract normalized counts matrix
  normalized_counts <- counts(dds_obj, normalized=TRUE) %>%
    as_tibble(rownames='genes')

  dds_object <- DESeq2::estimateSizeFactors(dds_obj)

  top_ten <- labeled_results %>%
    #slice_min = select 10 rows with lowest padj value
    dplyr::slice_min(padj, n = num_genes) %>%
    dplyr::select(genes)

  #join top 10 selected gene id with normalized count, matching values by genes 
  top_ten <- left_join(top_ten,normalized_counts, by='genes') %>%
    #collapse sample names and norm counts into one column // changes wide data to long data// exclude first column (genes)
    gather(sample_names, norm_counts, -1)
  
  scatter <- top_ten %>%
    ggplot2::ggplot() +
    geom_jitter(aes(x=genes, y=log10(norm_counts), color=sample_names)) +
    theme_light() + 
    ggtitle('Log10 Normalized Counts for Top Ten DE Genes') +
    #change angle of xlables to be readable
    theme(axis.text.x = element_text(angle = 90)) %>%
    return()
}


#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
  volcano <- labeled_results %>%
    ggplot2::ggplot(aes(x=log2FoldChange,y=-log10(padj),color=volc_plot_status)) +
    geom_point() +
    theme_light()+
    geom_hline(yintercept = 0.1, linetype = 'dashed') +
    ggtitle('Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)') %>%
    return()


}

#' Function to run fgsea on DESeq2 results
#'
#' @param labeled_results (tibble): the labeled results from DESeq2
#' @param gmt (str): the path to the GMT file
#' @param min_size: the threshold for minimum size of the gene set
#' @param max_size: the threshold for maximum size of the gene set
#'
#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
run_gsea <- function(labeled_results, gmt, min_size, max_size) {
  
  return(NULL)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
  return(NULL)
}

