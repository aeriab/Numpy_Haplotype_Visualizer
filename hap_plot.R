#!/usr/bin/env Rscript

# Peter Laurin and Brendan Aeria
# Modified Jan 2026 for Batch Processing Support
# General use haplotype plotting script

library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(magrittr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(argparse, warn.conflicts = FALSE)

# Detect support packages quietly
has_reticulate <- suppressMessages(require("reticulate", quietly = TRUE))
has_rcppcnpy <- suppressMessages(require("RcppCNPy", quietly = TRUE))

get_args <- function() {
  parser <- ArgumentParser(description = "General use haplotype plotting script")
  parser$add_argument("-f", "--file", type="character", required=TRUE, metavar="PATH",
                      help="input file name (should be .tsv, .csv, or .npy).")
  parser$add_argument("-o", "--out", type="character", default="hap_plot.png", metavar="PATH",
                      help="output file name or directory (if index is -2)")
  
  parser$add_argument("--sort_method", type="character", default="frequency",
                        help="sorting method: 'none', 'frequency', or 'distance'")
  
  parser$add_argument("--window", type="double", nargs=2, default=c(0, 0), metavar=c("START", "END"),
                      help="window region to plot / cluster haplotypes by")
  
  parser$add_argument("-i", "--image_index", type="integer", default=-1,
                      help="Image index. Set to -2 to process entire batch into a folder.")
  
  parser$add_argument("--annotate", action="store_true", default=FALSE,
                      help="annotate WINDOW region in plot")
  parser$add_argument("--palette", type="character", default="default",
                      help="color minor alleles: 'default' or 'site_type'")
  parser$add_argument("--expanded_region", type="double", nargs=2, default=c(0, 0), metavar=c("START", "END"))
  parser$add_argument("--nonsample_cols", type="character", nargs="*")
  parser$add_argument("--width", type="double", default=6)
  parser$add_argument("--height", type="double", default=4)
  parser$add_argument("--polarize_to_minor", action="store_true", default=FALSE)
  args <- parser$parse_args()
  return(args)
}

get_haplotype_clusters <- function(raw_haps){
  matrix_haps <- t(raw_haps)
  hap_dists <- dist(matrix_haps, method="manhattan")
  hap_clust <- hclust(hap_dists, method="average")
  raw_cluster <- cutree(hap_clust, h = 0)
  cluster_counts <- table(raw_cluster) %>% as_tibble() %>%
                    arrange(desc(n)) %>%
                    mutate(raw_cluster = as.integer(raw_cluster), ranked_hap = row_number())

  cluster_map <- tibble(sample = names(raw_cluster), raw_cluster) %>%
    left_join(cluster_counts,by='raw_cluster') %>%
    arrange(ranked_hap)

  dist_mat <- as.matrix(hap_dists, labels=TRUE)
  rep_hap <- cluster_map %>% filter(ranked_hap == 1) %>% pull(sample)
  rep_hap <- rep_hap[1]

  dist_map <- cluster_map %>% mutate(dist_to_rep = dist_mat[rep_hap,sample]) %>%
              group_by(ranked_hap) %>%
              summarize(ranked_hap=ranked_hap[1],dist_to_rep=dist_to_rep[1]) %>%
              mutate(ranked_hap_dist = rank(dist_to_rep,ties.method="first"))

  final_cluster_map <- cluster_map %>% left_join(dist_map, by="ranked_hap") %>%
                       rename(ranked_hap_freq = ranked_hap)
  return(final_cluster_map)
}

plot_haplotype <- function(l,r,hap,sort_method="frequency",nonsample_cols=NA,annotation=F,palette='default', polarize=F){
  sample_cols <- colnames(hap)[!colnames(hap) %in% nonsample_cols]
  if(polarize){
    afs <- rowMeans(hap %>% select(all_of(sample_cols)), na.rm=TRUE)
    to_polarize <- afs > 0.5
    hap[to_polarize, sample_cols] <- 1 - hap[to_polarize, sample_cols]
  }
  hap_df <- hap %>% mutate(site_index = rank(site_pos)) %>%
    pivot_longer(cols=-any_of(c(nonsample_cols, "site_index")), names_to = "sample", values_to = "base")
  
  hap_df$sample <- factor(hap_df$sample, levels = sample_cols)

  if (sort_method != "none") {
      raw_hap <- hap %>% filter(site_pos >= l & site_pos <= r) %>%
                         select(-any_of(nonsample_cols)) %>%
                         select(where(~ !all(is.na(.))))
      hap_clusters <- get_haplotype_clusters(raw_hap)
      if(sort_method == "frequency"){
        new_order <- hap_clusters %>% arrange(ranked_hap_freq) %>% pull(sample)
      } else if(sort_method == "distance"){
        new_order <- hap_clusters %>% arrange(ranked_hap_dist) %>% pull(sample)
      }
      hap_df$sample <- factor(hap_df$sample, levels = new_order)
  }
  
  hap_df$haplo_indices <- as.integer(hap_df$sample)
  n_samples <- length(unique(hap_df$sample))
  y_breaks <- seq(max(-n_samples,-10),-n_samples,by=-10)

  n_sites <- max(hap_df$site_index)
  n_breaks <- if(n_sites < 10) 1 else if(n_sites < 50) 5 else 10
  break_indices <- seq.int(1, max(hap_df$site_index), length.out=n_breaks) %>% round()
  x_breaks <- unique(hap_df$site_index)[break_indices]
  x_labels <- unique(hap_df$site_pos)[break_indices]
  pos_label <- "genomic position"
  if(max(x_labels) > 1e6){
    x_labels <- round(x_labels / 1e6,3); pos_label <- paste0(pos_label, " (Mb)")
  } else if (max(x_labels) > 1e3){
    x_labels <- round(x_labels / 1e3,3); pos_label <- paste0(pos_label, " (kb)")
  }

  if(palette == 'site_type'){
    hap_df <- hap_df %>% mutate(base = case_when(
      base == 0 ~ "0", is.na(base) ~ NA,
      site_type == 'syn' & base == 1 ~ "1",
      site_type == 'nonsyn' & base == 1 ~ "2",
      site_type %in% c("non_coding","NC") & base == 1 ~ "3",
      .default = "1"
    ))
  }

  plt <- ggplot(hap_df,aes(x=site_index,y=desc(haplo_indices), fill=as.character(base))) +
      geom_tile(show.legend = F) +
      theme(axis.title.y=element_blank(),panel.background=element_blank()) +
      scale_x_continuous(pos_label,breaks=x_breaks, labels=x_labels,expand=c(0.01,0)) +
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
  sort_mth <- args$sort_method
  file_ext <- tools::file_ext(args$file)
  nonsample_cols <- c("site_pos", "site_type", "contig", "gene_id")

  if (file_ext %in% c("tsv", "txt", "csv")) {
    if ((args$window[1] == 0) & (args$window[2] == 0)) stop("Error: --window required for TSV")
    delim_char <- if(file_ext == "csv") "," else "\t"
    hap <- read_delim(args$file, delim=delim_char, col_types = cols())
    if (!is.null(args$nonsample_cols)) nonsample_cols <- unique(c(nonsample_cols, args$nonsample_cols))
    l <- args$window[1]; r <- args$window[2]
    if(!((args$expanded_region[1] == 0) & (args$expanded_region[2] == 0))){
      hap <- hap %>% filter(site_pos >= args$expanded_region[1] & site_pos <= args$expanded_region[2])
    } else {
      hap <- hap %>% filter(site_pos >= l & site_pos <= r)
    }
    p <- plot_haplotype(l, r, hap, sort_method=sort_mth, nonsample_cols=nonsample_cols, annotation=args$annotate, palette=args$palette, polarize=args$polarize_to_minor)
    ggsave(args$out, plot=p, width=args$width, height=args$height)

  } else if (file_ext == "npy") {
    if(!has_reticulate && !has_rcppcnpy) stop("Error: Install 'RcppCNPy' or 'reticulate'.")
    
    if(has_reticulate){
      np <- import("numpy")
      full_array <- np$load(args$file)
    } else {
      full_array <- RcppCNPy::npyLoad(args$file)
    }

    # BATCH MODE: process all images into a folder
    if (args$image_index == -2) {
      if (!dir.exists(args$out)) dir.create(args$out, recursive = TRUE)
      batch_size <- dim(full_array)[1]
      for (idx in 1:batch_size) {
        raw_matrix <- full_array[idx, , ]
        hap <- as.data.frame(t(raw_matrix))
        colnames(hap) <- paste0("ind_", 1:ncol(hap))
        hap$site_pos <- 1:nrow(hap)
        hap <- hap %>% select(site_pos, everything())
        p <- plot_haplotype(1, nrow(hap), hap, sort_method=sort_mth, nonsample_cols=nonsample_cols, annotation=args$annotate, palette=args$palette, polarize=args$polarize_to_minor)
        ggsave(file.path(args$out, paste0("hap_idx_", idx, ".png")), plot=p, width=args$width, height=args$height)
      }
    } else {
      # SINGLE IMAGE MODE
      idx <- if(args$image_index == -1) 1 else args$image_index
      raw_matrix <- if(length(dim(full_array)) == 3) full_array[idx, , ] else full_array
      hap <- as.data.frame(t(raw_matrix))
      colnames(hap) <- paste0("ind_", 1:ncol(hap))
      hap$site_pos <- 1:nrow(hap)
      hap <- hap %>% select(site_pos, everything())
      p <- plot_haplotype(1, nrow(hap), hap, sort_method=sort_mth, nonsample_cols=nonsample_cols, annotation=args$annotate, palette=args$palette, polarize=args$polarize_to_minor)
      ggsave(args$out, plot=p, width=args$width, height=args$height)
    }
  }
}

main()
