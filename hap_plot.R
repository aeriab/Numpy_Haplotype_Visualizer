#!/usr/bin/env Rscript

# Peter Laurin and Brendan Aeria
# Nov 18, 2025 (Modified Nov 20, 2025 for numpy Support)
# General use haplotype plotting script

# -----------------------------------------------------------------------------
# Script Description:
#   Generates a haplotype image from either genomic TSV data or .npy image batches.
#   Auto-detects file format.
#
# Usage (TSV or CSV) example:
#   Rscript hap_plot.R --file data.tsv --window 10300000 10400000 --sort_method 'dist_sort'
#
# Usage (NPY) example:
#   Rscript hap_plot.R --file batch_data.npy --image_index 5 --sort_method 'cluster'
#
# Inputs:
#   -f / --file : delimited haplotype file or .npy file (Batch x Ind x Sites or Ind x Sites)
#   --window    : (Required for TSV) Start and End coordinates.
#   --image_index: (Required for NPY) The index of the image/batch to plot.
# -----------------------------------------------------------------------------

library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(argparse, warn.conflicts = FALSE)

get_args <- function() {
  parser <- ArgumentParser(description = "General use haplotype plotting script")
  parser$add_argument("-f", "--file", type="character", required=TRUE, metavar="PATH",
                      help="input file name (should be .tsv, .csv, or .npy).")
  parser$add_argument("-o", "--out", type="character", default="hap_plot.png", metavar="PATH",
                      help="output file name; extension should be one of '.png' '.jpeg' '.pdf'")
  
  parser$add_argument("--sort_method", type="character", default="frequency",
                        help="sorting method: 'none', 'frequency', or 'distance'")
  
  # Window is strictly required for TSV, ignored for NPY
  # PJL -- in a previous version this was not required -- i.e. you could plot the 
  # entire hap file, by default. Can we revert to this version? 
  parser$add_argument("--window", type="double", nargs=2, default=c(0, 0), metavar=c("START", "END"),
                      help="window region to plot / cluster haplotypes by")
  
  # Image index is strictly required for NPY, ignored for TSV
  # PJL -- likewise, can we make this optional, defaulting to an assumption that
  # the user can pass in a 2D array of a single haplotype window, rather than a 
  # long array of different windows? I've started implementing this below. 
  parser$add_argument("-i", "--image_index", type="integer", default=-1,
                      help="Image batch index. Required for .npy")
  
  parser$add_argument("--annotate", action="store_true", default=FALSE,
                      help="annotate WINDOW region in plot")
  parser$add_argument("--palette", type="character", default="default",
                      help="how to color minor alleles: 'default' or 'site_type' ('site_type' column must be in input.)")
  parser$add_argument("--expanded_region", type="double", nargs=2, default=c(0, 0), metavar=c("START", "END"),
                      help="expanded region to plot (and not cluster by)")
  parser$add_argument("--nonsample_cols", type="character", nargs="*",
                      help="non-sample columns (space separated). Commonly used columns for our lab (e.g. contig, gene_id) are automatically included.")
  parser$add_argument("--width", type="double", default=6,
                      help="output figure width")
  parser$add_argument("--height", type="double", default=4,
                      help="output figure height")
  parser$add_argument("--polarize_to_minor", action="store_true", default=FALSE,
                      help="polarize minor alleles to 1, major alleles to 0")
  args <- parser$parse_args()
  return(args)
}

get_haplotype_clusters <- function(raw_haps){

  # get hamming (manhattan) distances between haplotypes, 
  # and do hierarchical clustering (average ~ UPGMA)

  matrix_haps <- t(raw_haps)
  
  hap_dists <- dist(matrix_haps, method="manhattan")
  hap_clust <- hclust(hap_dists, method="average")

  # get haplotype clusters that cluster perfectly (h=0)
  raw_cluster <- cutree(hap_clust, h = 0)
  cluster_counts <- table(raw_cluster) %>% as_tibble() %>%
                    arrange(desc(n)) %>%
                    mutate(raw_cluster = as.integer(raw_cluster),
                           ranked_hap = row_number())


  # map samples to clusters, ranked by frequency
  cluster_map <- tibble(sample = names(raw_cluster), raw_cluster) %>%
    left_join(cluster_counts,by='raw_cluster') %>%
    arrange(ranked_hap)


  # get clusters' distance to most frequent haplotype cluster
  dist_mat <- as.matrix(hap_dists, labels=TRUE)
  rep_hap <- cluster_map %>% filter(ranked_hap == 1) %>% pull(sample)
  rep_hap <- rep_hap[1]

  dist_map <- cluster_map %>% mutate(dist_to_rep = dist_mat[rep_hap,sample]) %>%
              group_by(ranked_hap) %>%
              summarize(ranked_hap=ranked_hap[1],dist_to_rep=dist_to_rep[1]) %>%
              mutate(ranked_hap_dist = rank(dist_to_rep,ties.method="first"))

  # join all together
  final_cluster_map <- cluster_map %>% left_join(dist_map, by="ranked_hap") %>%
                       rename(ranked_hap_freq = ranked_hap)

  return(final_cluster_map)
}

#

plot_haplotype <- function(l,r,hap,sort_method="frequency",nonsample_cols=NA,annotation=F,palette='default', polarize=F){
  
  # 1. Identify Sample Columns (preserve original order from file)
  #    This ensures that 'ind_2' comes after 'ind_1', not 'ind_10'
  sample_cols <- colnames(hap)[!colnames(hap) %in% nonsample_cols]
  
  # 2. Pivot Data
  if(polarize){
    afs <- rowMeans(hap %>% select(all_of(sample_cols)), na.rm=TRUE)
    to_polarize <- afs > 0.5
    hap[to_polarize, sample_cols] <- 1 - hap[to_polarize, sample_cols]
  }
  hap_df <- hap %>% mutate(site_index = rank(site_pos)) %>%
    pivot_longer(cols=-any_of(c(nonsample_cols, "site_index")), names_to = "sample", values_to = "base")
  
  # 3. Default: Set Factor Levels to Original File Order ("none")
  hap_df$sample <- factor(hap_df$sample, levels = sample_cols)

  # 4. Conditional: Update Factor Levels if Clustering is requested
  if (sort_method != "none") {
      
      raw_hap <- hap %>% filter(site_pos >= l & site_pos <= r) %>% 
                         select(-any_of(nonsample_cols)) %>% 
                         select(where(~ !all(is.na(.)))) # filter out all-NA columns
                                
      hap_clusters <- get_haplotype_clusters(raw_hap)
      
      # We don't strictly need to join the whole table, just get the order
      if(sort_method == "frequency"){
        new_order <- hap_clusters %>% arrange(ranked_hap_freq) %>% pull(sample)
      } else if(sort_method == "distance"){
        new_order <- hap_clusters %>% arrange(ranked_hap_dist) %>% pull(sample)
      }
      
      # Re-level the factor to the new sorted order
      hap_df$sample <- factor(hap_df$sample, levels = new_order)
  }
  
  # 5. Generate Y-axis indices based on the Factor Integer Value
  #    (This guarantees the plot follows the factor levels we set above)
  hap_df$haplo_indices <- as.integer(hap_df$sample)

  # get axes labels
  n_samples <- length(unique(hap_df$sample))
  # ... (rest of function below remains exactly the same)
  y_breaks <- seq(max(-n_samples,-10),-n_samples,by=-10)

  n_sites <- max(hap_df$site_index)
  if(n_sites < 10){
    n_breaks <- 1
  } else if (n_sites < 50) {
    n_breaks <- 5
  } else {
    n_breaks <- 10
  }
  break_indices <- seq.int(1, max(hap_df$site_index), length.out=n_breaks) %>% round()
  x_breaks <- unique(hap_df$site_index)[break_indices]
  x_labels <- unique(hap_df$site_pos)[break_indices]
  pos_label <- "genomic position"
  if(max(x_labels) > 1e6){
    x_labels <- round(x_labels / 1e6,3)
    pos_label <- paste0(pos_label, " (Mb)")
  } else if (max(x_labels) > 1e3){
    x_labels <- round(x_labels / 1e3,3)
    pos_label <- paste0(pos_label, " (kb)")
  }

  if(palette == 'site_type'){
    hap_df <- hap_df %>% mutate(base = case_when(
      base == 0 ~ "0",
      is.na(base) ~ NA, 
      site_type == 'syn' & base == 1 ~ "1",
      site_type == 'nonsyn' & base == 1 ~ "2",
      site_type %in% c("non_coding","NC") & base == 1 ~ "3",
      .default = "1"
    ))
  }

  # make plot
  plt <- ggplot(hap_df,aes(x=site_index,y=desc(haplo_indices), fill=as.character(base)))
  plt <- plt + geom_tile(show.legend = F)
  plt <- plt + theme(axis.title.y=element_blank(),panel.background=element_blank()) +
      scale_x_continuous(pos_label,breaks=x_breaks,
                         labels=x_labels,expand=c(0.01,0)) +
      scale_y_continuous(breaks=y_breaks,labels=abs, expand=c(0.01,0))
    
  if(palette == 'default'){
    plt <- plt + scale_fill_manual(values=c("0"="grey87","1"="steelblue2"),na.value="white")
  } else if (palette == 'site_type'){
    plt <- plt + scale_fill_manual(values=c("0"="grey87","1"="steelblue2","2"="firebrick3","3"="orange2"),na.value="white")
  }

  if(annotation){
    rect_l <- which.max(unique(hap_df$site_pos)>=l)
    rect_r <- which.max(unique(hap_df$site_pos)>r)-1
    plt <- plt + annotate("rect", xmin=rect_l, xmax=rect_r, ymin=0, ymax=n_samples %/% 100 + 1, fill="red")
  }

  return(plt)
}




main <- function(){

  args = get_args()
  
  sort_mth <- args$sort_method # "none", "frequency", "distance"
  
  # Detect file extension
  file_ext <- tools::file_ext(args$file)
  
  hap <- NULL
  l <- 0
  r <- 0
  
  # handle nonsample_cols
  nonsample_cols <- c("site_pos", "site_type", "contig", "gene_id")

  # ---------------------
  # MODE 1: delimited file
  # ---------------------
  if (file_ext %in% c("tsv", "txt", "csv")) {
    
    # Enforce Window Argument
    if ((args$window[1] == 0) & (args$window[2] == 0)) {
      stop("Error: For .tsv files, you must specify a --window (e.g., --window 10000 20000)")
    }
    
    if(file_ext == "tsv" || file_ext == "txt"){
      delim_char <- "\t"
    } else if (file_ext == "csv"){
      delim_char <- ","
    }
    hap <- read_delim(args$file, delim=delim_char, col_types = cols())
    
    # Handle nonsample columns
    if (!is.null(args$nonsample_cols)) {
      nonsample_cols <- unique(c(nonsample_cols, args$nonsample_cols))
    }
    
    # Window Logic
    l <- args$window[1]
    r <- args$window[2]
    
    expanded_region <- args$expanded_region
    if(!((expanded_region[1] == 0) & (expanded_region[2] == 0))){
      hap <- hap %>% filter(site_pos >= expanded_region[1] & site_pos <= expanded_region[2])
    } else {
      hap <- hap %>% filter(site_pos >= l & site_pos <= r)
    }

  # ---------------------
  # MODE 2: .npy File
  # ---------------------
  } else if (file_ext == "npy") {
    
    # PJL -- move here, don't require average user to install reticulate / Python
    # rccppcnpy is lighter weight, but only handles 2d arrays
    # assume user can pass in a 2d numpy array for single image
    if(!require("RcppCNPy") && !require("reticulate")){ 
      stop("Error: To read .npy files, please install either the 'RcppCNPy' or 'reticulate' package.")
    }

    if(require("reticulate")){
      np <- import("numpy")
      full_array <- np$load(args$file)
    } else if (require("RcppCNPy") & args$image_index == -1){
      full_array <- RcppCNPy::npyLoad(args$file, type="integer") # needs to be 64-bit
    }

    if(args$image_index == -1){
      raw_matrix <- full_array
    } else if(args$image_index <= dim(full_array)[1] || args$image_index >= 1) {
      raw_matrix <- full_array[args$image_index, , ]
    } else {
      stop(paste("Image index", args$image_index, "is out of bounds. Max batch size is", dim(full_array)[1]))
    }

    
    # Transpose to match pipeline: Rows=Sites, Cols=Individuals
    hap_matrix <- t(raw_matrix)
    hap <- as.data.frame(hap_matrix)
    
    # Synthesize columns
    # PJL -- I like this default columns approach -- can we implement in TSV mode too?
    colnames(hap) <- paste0("ind_", 1:ncol(hap))
    hap$site_pos <- 1:nrow(hap)
    hap <- hap %>% select(site_pos, everything())
    
    # Set window to full image
    l <- 1
    r <- nrow(hap)
    
  } else {
    stop(paste("Unsupported file extension:", file_ext, ". Please use .tsv or .npy"))
  }

  # ---------------------
  # PLOTTING
  # ---------------------
  p <- plot_haplotype(l, r, hap, sort_method=sort_mth, nonsample_cols=nonsample_cols, annotation=args$annotate, palette=args$palette, polarize=args$polarize_to_minor)
  ggsave(args$out, plot=p, width=args$width, height=args$height)
}

main()
