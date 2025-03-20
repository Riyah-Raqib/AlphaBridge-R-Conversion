# ==============================
# Load Required Libraries
# ==============================
library(R6)
library(Biostrings)
library(jsonlite)

protein_seq_alingment <- R6Class("protein_seq_alingment",
  public = list(
    query_seq = NULL,
    ref_seq = NULL,
    
    initialise = function(query_seq, ref_seq) {
      self$query_seq <- query_seq
      self$ref_seq <- ref_seq
    },
    
    get_alignment = function() {
      ref_seq <- self$query_seq
      query_seq <- self$ref_seq
      
      # Create aligner equivalent using Biostrings::pairwiseAlignment with BLASTP-like scoring
      # Convert sequences to AAString objects
      aa_ref <- AAString(ref_seq)
      aa_query <- AAString(query_seq)
      
      # Set alignment parameters: global alignment with BLOSUM62 matrix.
      # Gap penalties are set similar to typical BLASTP defaults.
      alignment <- pairwiseAlignment(aa_ref, aa_query,
                                     substitutionMatrix = BLOSUM62,
                                     gapOpening = 11,
                                     gapExtension = 1,
                                     type = "global")
      return(alignment)
    }
  )
)

compare_protein_seq <- R6Class("compare_protein_seq",
  public = list(
    structure_list = NULL,
    fasta_list = NULL,
    
    initialise = function(structure_list, fasta_list) {
      self$structure_list <- structure_list
      self$fasta_list <- fasta_list
    },
    
    extract_chain_dict = function() {
      structure_list <- self$structure_list
      fasta_list <- self$fasta_list
      
      chain_dict <- list()
      
      # Using lists to mimic defaultdict behavior
      struct_dict <- list()
      fasta_dict <- list()
      chain_dict <- list()
      
      n <- min(length(structure_list), length(fasta_list))
      for(i in seq_len(n)) {
        struct <- structure_list[[i]]
        fasta <- fasta_list[[i]]
        
        # Unpack struct and fasta assuming each is a two-element vector/list
        struct_chain <- struct[[1]]
        struct_seq <- struct[[2]]
        fasta_chain <- fasta[[1]]
        fasta_seq <- fasta[[2]]
        
        # If not already in struct_dict, initialise as empty vector
        if(is.null(struct_dict[[as.character(struct_seq)]])) {
          struct_dict[[as.character(struct_seq)]] <- c()
        }
        if(is.null(fasta_dict[[as.character(fasta_seq)]])) {
          fasta_dict[[as.character(fasta_seq)]] <- c()
        }
        
        struct_dict[[as.character(struct_seq)]] <- c(struct_dict[[as.character(struct_seq)]], struct_chain)
        fasta_dict[[as.character(fasta_seq)]] <- c(fasta_dict[[as.character(fasta_seq)]], fasta_chain)
      }
      
      key_seq_union <- union(names(struct_dict), names(fasta_dict))
      
      for(key_seq in key_seq_union) {
        vec_struct <- struct_dict[[key_seq]]
        vec_fasta <- fasta_dict[[key_seq]]
        n2 <- min(length(vec_struct), length(vec_fasta))
        for(j in seq_len(n2)) {
          chain_dict[[ vec_fasta[j] ]] <- vec_struct[j]
        }
      }
      
      return(chain_dict)
    }
  )
)

msa_folder <- R6Class("msa_folder",
  public = list(
    feature_folder = NULL,
    fasta_list = NULL,
    
    initialise = function(feature_folder, fasta_list) {
      self$feature_folder <- feature_folder
      self$fasta_list <- fasta_list
    },
    
    extract_chain_id_map = function() {
      feature_folder <- self$feature_folder
      chain_id_map_path <- file.path(feature_folder, "msas", "chain_id_map.json")
      
      json_file <- readLines(chain_id_map_path, warn = FALSE)
      chain_id_map <- fromJSON(paste(json_file, collapse=""))
      
      return(chain_id_map)
    },
    
    extract_msa_folder = function() {
      feature_folder <- self$feature_folder
      chain_id_map <- self$extract_chain_id_map()
      
      parent_msa_folder <- file.path(feature_folder, "msas")
      
      dedup_msa_folder <- list()
      
      for(record in names(chain_id_map)) {
        sequence <- chain_id_map[[record]]$sequence
        msa_folder_path <- file.path(parent_msa_folder, record)
        if(file.exists(msa_folder_path)) {
          dedup_msa_folder[[ as.character(sequence) ]] <- record
        }
      }
      
      return(dedup_msa_folder)
    },
    
    map_description2msa_folder = function() {
      dedup_msa_folder <- self$extract_msa_folder()
      fasta_list <- self$fasta_list
      
      fasta_msa_dict <- list()
      
      n <- length(fasta_list)
      for(i in seq_len(n)) {
        fasta_entry <- fasta_list[[i]]
        fasta_name <- fasta_entry[[1]]
        fasta_seq <- fasta_entry[[2]]
        
        if(is.null(fasta_msa_dict[[fasta_name]])) {
          fasta_msa_dict[[fasta_name]] <- dedup_msa_folder[[ as.character(fasta_seq) ]]
        }
      }
      
      return(fasta_msa_dict)
    }
  )
)
