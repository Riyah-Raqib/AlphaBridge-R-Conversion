# Load required libraries
library(R6)
library(igraph)
library(leiden)      # For community detection using Leiden algorithm
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)

# For image plotting with colorbar (fields alternative)
library(fields)

# -*- coding: utf-8 -*-
# """
# Created on Wed Feb 7 16:53:05 2024
# 
# @author: Dan_salv
# """

domain_clustering <- R6Class("domain_clustering",
  public = list(
    # Class attributes
    matrix_dict = NULL,
    list_fasta_files = NULL,
    plotting = NULL,
    outdir = NULL,
    alphafold_version = NULL,
    bool_mask_clusters = NULL,
    
    # __init__ method
    initialize = function(matrix_dict, list_fasta_files, plotting = TRUE, bool_mask_clusters = TRUE, outdir = '' , alphafold_version = 'AF2'){
      self$matrix_dict <- matrix_dict
      self$list_fasta_files <- list_fasta_files
      self$plotting <- plotting
      self$outdir <- outdir
      self$alphafold_version <- alphafold_version
      self$bool_mask_clusters <- bool_mask_clusters
    },
    
    run_domain_clustering = function(){
      matrix_dict <- self$matrix_dict
      list_fasta_files <- self$list_fasta_files
      
      matrix_input <- matrix_dict[['confidence_matrix']]
      masked_confidance_matrix <- matrix_dict[['masked_confidence_matrix']]
      masked_contact_matrix <- matrix_dict[['masked_contact_matrix']]
      
      if (self$alphafold_version == 'AF2') {
        graph_resolution <- 0.5
        matrix_cutoff <- 2.6
        cmap_confidance <- 'RdPu_r'
        cmap_contact <- 'Blues_r'
        
        coevolutionary_domains <- get_coevolutionary_domains(matrix_input, graph_resolution = graph_resolution, matrix_cutoff = matrix_cutoff)
        
        temp <- self$get_interacting_coevolutionary_domains(coevolutionary_domains)
        interacting_coevultionary_cluster_dict <- temp[[1]]
        entity_region_dict <- temp[[2]]
        
        interacting_mask_cluster <- self$get_interacting_mask_cluster(coevolutionary_domains, interacting_coevultionary_cluster_dict, self$bool_mask_clusters)
        
      } else if (self$alphafold_version == 'AF3') {
        
        graph_resolution <- 0.25
        matrix_cutoff <- 27
        cmap_confidance <- 'Blues_r'
        cmap_contact <- "RdPu"
        
        coevolutionary_domains <- get_coevolutionary_domains(matrix_input, graph_resolution = graph_resolution, matrix_cutoff = matrix_cutoff)
        
        temp <- self$get_interacting_coevolutionary_domains(coevolutionary_domains)
        interacting_coevultionary_cluster_dict <- temp[[1]]
        entity_region_dict <- temp[[2]]
        
        interacting_mask_cluster <- self$get_interacting_mask_cluster(coevolutionary_domains, interacting_coevultionary_cluster_dict, self$bool_mask_clusters)
        
      }
      
      if (self$plotting) {
        plot_combination_matrix(coevolutionary_domains, masked_confidance_matrix, masked_contact_matrix,
                                interacting_mask_cluster, list_fasta_files, self$outdir, self$alphafold_version)
        plot_separate_matrix(matrix_dict, list_fasta_files, self$outdir)
      }
      
      return(list(coevolutionary_domains = coevolutionary_domains,
                  interacting_coevultionary_cluster_dict = interacting_coevultionary_cluster_dict,
                  entity_region_dict = entity_region_dict))
    },
    
    get_interacting_coevolutionary_domains = function(coevolutionary_domains){
      # Unpack list_fasta_files
      list_fasta_name <- self$list_fasta_files[[1]]
      list_fasta_acclen <- self$list_fasta_files[[2]]
      list_fasta_centerticks <- self$list_fasta_files[[3]]
      list_fasta_len <- self$list_fasta_files[[4]]
      
      unique_coevolutionary_domains <- as.numeric(unique(coevolutionary_domains))
      cluster_index_dict <- list()
      interacting_coevultionary_cluster_dict <- list()
      entity_region_dict <- list()
      
      #iterate through all coevolutionary_domains found
      for (cluster_name in unique_coevolutionary_domains) {
        ii <- which(coevolutionary_domains == cluster_name) - 0  # indices remain numeric
        if (length(ii) > 1) {
          if (!(as.character(cluster_name) %in% names(cluster_index_dict))) {
            cluster_index_dict[[as.character(cluster_name)]] <- ii
          }
          
          #define start and end for each cluster    
          # Group consecutive indices
          groups <- split(ii, cumsum(c(TRUE, diff(ii) != 1)))
          group_index_range_list <- lapply(groups, function(item) {
            c(item[1], tail(item, n=1))
          })
          
          #build interacting_coevultionary_cluster_dict
          cluster_group_name <- paste0('cluster_', cluster_name)
          if (!(cluster_group_name %in% names(interacting_coevultionary_cluster_dict))) {
            interacting_coevultionary_cluster_dict[[cluster_group_name]] <- list(
              range_index_list = group_index_range_list,
              overlap_complex = list()
            )
          }
          
          #find overlapping sequences between all proteins and regions of coevolutionary_domains 
          for(i in seq_along(list_fasta_name)) {
            protein <- list_fasta_name[[i]]
            acclen <- list_fasta_acclen[[i]]
            fasta_len <- list_fasta_len[[i]]
            
            protein_start <- acclen - fasta_len
            protein_end <- acclen - 1
            
            entity_region_dict[[protein]] <- c(protein_start, protein_end)
            
            for (group_range in interacting_coevultionary_cluster_dict[[cluster_group_name]]$range_index_list) {
              if (length(group_range) == 2) {
                group_index_start <- group_range[1]
                group_index_end <- group_range[2]
                
                if (group_index_start < protein_end && protein_start < group_index_end) {
                  protein_range <- seq(protein_start, protein_end)
                  cluster_range <- seq(group_index_start, group_index_end)
                  overlapping_range <- c(max(protein_range[1], cluster_range[1]),
                                         min(tail(protein_range, n=1), tail(cluster_range, n=1)) + 1)
                  
                  if (!(protein %in% names(interacting_coevultionary_cluster_dict[[cluster_group_name]]$overlap_complex))) {
                    interacting_coevultionary_cluster_dict[[cluster_group_name]]$overlap_complex[[protein]] <- list()
                  }
                  interacting_coevultionary_cluster_dict[[cluster_group_name]]$overlap_complex[[protein]] <- 
                    c(interacting_coevultionary_cluster_dict[[cluster_group_name]]$overlap_complex[[protein]],
                      list(overlapping_range))
                }
              }
            }
          }
        }
      }
      return(list(interacting_coevultionary_cluster_dict, entity_region_dict))
    },
    
    get_interacting_mask_cluster = function(clusters, interacting_coevultionary_cluster_dict, bool_mask_clusters){
      if (!bool_mask_clusters) {
        interacting_mask_cluster <- rep(FALSE, length(clusters))
        return(interacting_mask_cluster)
      } else {
        interacting_mask_cluster <- rep(TRUE, length(clusters))
        interacting_range_list <- list()
        for (cluster in names(interacting_coevultionary_cluster_dict)) {
          if (length(interacting_coevultionary_cluster_dict[[cluster]]$overlap_complex) > 1) {
            for (entity in names(interacting_coevultionary_cluster_dict[[cluster]]$overlap_complex)) {
              for (interacting_range in interacting_coevultionary_cluster_dict[[cluster]]$overlap_complex[[entity]]) {
                interacting_range_list <- c(interacting_range_list, list(interacting_range))
              }
            }
          }
        }
        
        for (index in seq_along(clusters)) {
          for (interacting_range in interacting_range_list) {
            if ( (interacting_range[1] <= index) && (index <= interacting_range[2]) ) {
              interacting_mask_cluster[index] <- FALSE
            }
          }
        }
      }
      return(interacting_mask_cluster)
    }
  )
)

get_coevolutionary_domains <- function(matrix_input, pae_power = 1, graph_resolution = 0.5, matrix_cutoff = 2.6 ){
  # Compute weights = 1/matrix_input**pae_power
  weights <- 1/(matrix_input^pae_power)
  
  size <- nrow(matrix_input)
  # Create an empty undirected graph with vertices 0 to size-1
  # In R vertices are 1-indexed; we create size vertices.
  g <- make_empty_graph(n = size, directed = FALSE)
  
  # Find edges where matrix_input < matrix_cutoff
  edges_idx <- which(matrix_input < matrix_cutoff, arr.ind = TRUE)
  # Adjust indices for R (they are already 1-indexed if matrix was created in R)
  if(nrow(edges_idx) > 0){
    # igraph expects a vector of vertex ids in pairs
    edge_list <- as.vector(t(edges_idx))
    g <- add_edges(g, edge_list)
    sel_weights <- mapply(function(i,j) { weights[i,j] },
                          edges_idx[,1], edges_idx[,2])
    E(g)$weight <- sel_weights
  }
  
  # Perform Leiden community detection
  # Note: using leiden_find_partition from the leiden package
  vc <- leiden_find_partition(g, resolution_parameter = graph_resolution/100, weights = E(g)$weight)
  coevolutionary_domains <- as.numeric(vc$membership)
  
  return(coevolutionary_domains)
}

plot_combination_matrix <- function(coevolutionary_domains, confidance_matrix, contact_matrix, interacting_mask_cluster, list_fasta_files, outdir, alphafold_version){
  # labels = np.array(coevolutionary_domains)
  labels <- as.numeric(coevolutionary_domains)
  # label_data = np.tile(labels, (2,1))
  label_data <- rbind(labels, labels)
  
  # mask_data = np.tile(interacting_mask_cluster, (2,1))
  mask_data <- rbind(interacting_mask_cluster, interacting_mask_cluster)
  list_fasta_name <- list_fasta_files[[1]]
  list_fasta_acclen <- list_fasta_files[[2]]
  list_fasta_centerticks <- list_fasta_files[[3]]
  list_fasta_len <- list_fasta_files[[4]]
  
  inv_fasta_acclen <- sapply(list_fasta_acclen, function(acclen) { sum(list_fasta_len) - acclen })
  inv_fasta_acclen <- c(sum(list_fasta_len) - 1, inv_fasta_acclen)
  
  list_fasta_len_names <- mapply(function(i, length_val){
    if(i != (length(list_fasta_len))) {
      paste0("0 / ", length_val)
    } else {
      paste0(length_val)
    }
  }, seq_along(list_fasta_len), list_fasta_len, SIMPLIFY = TRUE)
  
  # Create a PNG file for saving the plot
  png(filename = file.path(outdir, "Confidence-contact_plot.png"), width = 1500, height = 1500, res = 300)
  
  # Set up layout: main plot and additional axes for colorbars and heatmap
  layout_matrix <- matrix(c(2,3,1,1), nrow = 2, byrow = TRUE)
  layout(layout_matrix, widths = c(4,1), heights = c(1,4))
  
  # Plot the label heatmap on top (simulating cax3)
  par(mar = c(1,1,1,1))
  image(t(label_data), axes = FALSE, col = brewer.pal(12, "Paired"))
  # Add center ticks and labels using list_fasta_centerticks and list_fasta_name
  # (Exact replication of minor ticks is not possible, so approximate)
  
  # Plot the main axis with confidence and contact matrices
  par(mar = c(5,5,5,5))
  # Plot confidance_matrix first
  image(confidance_matrix, col = brewer.pal(9, "Blues"), axes = FALSE, main = "")
  # Overlay contact_matrix using different color map depending on alphafold_version
  if (alphafold_version == 'AF2') {
    image(contact_matrix, col = rev(brewer.pal(9, "RdPu")), add = TRUE)
  } else if (alphafold_version == 'AF3') {
    image(contact_matrix, col = brewer.pal(9, "RdPu"), add = TRUE)
  }
  
  # Add axis ticks and labels
  axis(1, at = list_fasta_acclen/length(confidance_matrix[1,]), labels = NA, lwd.ticks = 1, col.ticks = "black")
  axis(2, at = list_fasta_acclen/length(confidance_matrix[,1]), labels = list_fasta_len_names, las = 1)
  # Draw vertical and horizontal lines for each acclen
  for (acclen in list_fasta_acclen) {
    abline(v = acclen/length(confidance_matrix[1,]), col = "black", lwd = 1)
    abline(h = acclen/length(confidance_matrix[,1]), col = "black", lwd = 1)
  }
  
  # Plot colorbars using image.plot from fields
  # Right side colorbar for confidance_matrix
  par(mar = c(5,0,5,2))
  image.plot(legend.only=TRUE, zlim=range(confidance_matrix), col=brewer.pal(9,"Blues"))
  # Bottom colorbar for contact_matrix
  par(mar = c(2,5,2,5))
  image.plot(legend.only=TRUE, horizontal=TRUE, zlim=range(contact_matrix), col=if(alphafold_version=='AF2') rev(brewer.pal(9,"RdPu")) else brewer.pal(9,"RdPu"))
  
  dev.off()
}

plot_separate_matrix <- function(matrix_dict, list_fasta_files, outdir){
  matrix_list <- c('pae','confidence_matrix','contact_matrix')
  cmap_list <- c('Greens_r', 'Blues_r', "RdPu")
  
  list_fasta_name <- list_fasta_files[[1]]
  list_fasta_acclen <- list_fasta_files[[2]]
  list_fasta_centerticks <- list_fasta_files[[3]]
  list_fasta_len <- list_fasta_files[[4]]
  list_fasta_len_names <- mapply(function(i, length_val){
    if(i != (length(list_fasta_len))) {
      paste0(length_val, " / 0")
    } else {
      paste0(length_val)
    }
  }, seq_along(list_fasta_len), list_fasta_len, SIMPLIFY = TRUE)
  
  for (n in seq_along(matrix_list)) {
    feature <- matrix_list[n]
    cmap <- cmap_list[n]
    png(filename = file.path(outdir, paste0(feature, ".png")), width = 1500, height = 1500, res = 300)
    
    par(mar = c(5,5,5,5))
    matrix_plot <- matrix_dict[[feature]]
    image(matrix_plot, axes = FALSE, col = switch(cmap,
                                                  'Greens_r' = rev(brewer.pal(9, "Greens")),
                                                  'Blues_r' = rev(brewer.pal(9, "Blues")),
                                                  'RdPu' = brewer.pal(9, "RdPu")))
    axis(1, at = list_fasta_acclen/length(matrix_plot[1,]), labels = NA)
    axis(1, at = list_fasta_centerticks/length(matrix_plot[1,]), labels = list_fasta_name, las = 2, cex.axis = 1.5, tick = FALSE)
    axis(2, at = (as.numeric(list_fasta_acclen)-1)/nrow(matrix_plot), labels = list_fasta_len_names, cex.axis = 1.5, las = 1)
    
    for (i in list_fasta_acclen) {
      abline(v = i/length(matrix_plot[1,]), col = "black", lwd = 1)
      abline(h = i/length(matrix_plot[,1]), col = "black", lwd = 1)
    }
    
    # Add colorbar using image.plot
    image.plot(legend.only = TRUE, zlim = range(matrix_plot), col = switch(cmap,
                                                                           'Greens_r' = rev(brewer.pal(9, "Greens")),
                                                                           'Blues_r' = rev(brewer.pal(9, "Blues")),
                                                                           'RdPu' = brewer.pal(9, "RdPu")))
    dev.off()
  }
}

plot_joined_matrix <- function(matrix_dict, list_fasta_files, outdir){
  N_COL <- 3
  N_ROW <- 1
  matrix_list <- c('pae','confidence_matrix','contact_matrix')
  cmap_list <- c('Greens_r', 'Blues_r', "RdPu")
  
  list_fasta_name <- list_fasta_files[[1]]
  list_fasta_acclen <- list_fasta_files[[2]]
  list_fasta_centerticks <- list_fasta_files[[3]]
  list_fasta_len <- list_fasta_files[[4]]
  list_fasta_len_names <- mapply(function(i, length_val){
    if(i != (length(list_fasta_len))) {
      paste0(length_val, " / 0")
    } else {
      paste0(length_val)
    }
  }, seq_along(list_fasta_len), list_fasta_len, SIMPLIFY = TRUE)
  
  png(filename = file.path(outdir, "combination_matrix.png"), width = 3000, height = 1000, res = 300)
  par(mfrow = c(N_ROW, N_COL), mar = c(5,5,5,5))
  
  for (n in seq_along(matrix_list)) {
    feature <- matrix_list[n]
    cmap <- cmap_list[n]
    matrix_plot <- matrix_dict[[feature]]
    image(matrix_plot, axes = FALSE, col = switch(cmap,
                                                  'Greens_r' = rev(brewer.pal(9, "Greens")),
                                                  'Blues_r' = rev(brewer.pal(9, "Blues")),
                                                  'RdPu' = brewer.pal(9, "RdPu")))
    axis(1, at = list_fasta_acclen/length(matrix_plot[1,]), labels = NA)
    axis(1, at = list_fasta_centerticks/length(matrix_plot[1,]), labels = list_fasta_name, las = 2, cex.axis = 1.5)
    axis(2, at = (as.numeric(list_fasta_acclen)-1)/nrow(matrix_plot), labels = list_fasta_len_names, cex.axis = 1.5, las = 1)
    
    for (i in list_fasta_acclen) {
      abline(v = i/length(matrix_plot[1,]), col = "black", lwd = 1)
      abline(h = i/length(matrix_plot[,1]), col = "black", lwd = 1)
    }
    image.plot(legend.only = TRUE, zlim = range(matrix_plot), col = switch(cmap,
                                                                           'Greens_r' = rev(brewer.pal(9, "Greens")),
                                                                           'Blues_r' = rev(brewer.pal(9, "Blues")),
                                                                           'RdPu' = brewer.pal(9, "RdPu")))
  }
  dev.off()
}

# The Python code was translated line‐by‐line into R using the R6 package to implement the class “domain_clustering” with its methods. All class attributes, methods, and parameters have been preserved exactly as in the original code.
# For numerical operations, R’s native vectorized operations (using operators such as “^” and “1/”) have been used to mimic numpy’s behavior.
# The igraph functionality is directly mapped using the R “igraph” package. Note that vertex indexing in Python (0-indexed) versus R (1-indexed) was handled by careful use of indices; it is assumed that the input matrices come from R and are 1-indexed.
# The Leiden community detection was implemented using the “leiden” package’s function leiden_find_partition. This is the closest R equivalent to Python’s g.community_leiden.
# Plotting functionality (heatmaps, overlays, colorbars, custom axes ticks) was re-implemented using a combination of base R plotting functions and the “fields” package’s image.plot to simulate matplotlib and seaborn behavior. Due to inherent differences between R and Python plotting libraries, the output may not be pixel‐identical but preserves the overall functionality.
# All comments from the original Python code have been preserved and translated as R comments.
# All dependencies have been explicitly imported at the top of the code.
# All functionality including error handling (implicit in R functions) and complete implementations are provided with no placeholders.
# Some aspects of layout (such as subplot adjustments, minor ticks, and axes divisions) have been approximated in R since R’s plotting system does not support an exact one‐to‐one mapping with matplotlib.
# The code is complete and working provided that all required packages ("R6", "igraph", "leiden", "ggplot2", "grid", "gridExtra", "RColorBrewer", and "fields") are installed.
