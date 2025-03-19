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
    
    initialize = function(query_seq, ref_seq) {
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
    
    initialize = function(structure_list, fasta_list) {
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
        
        # If not already in struct_dict, initialize as empty vector
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
    
    initialize = function(feature_folder, fasta_list) {
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

# Language-specific features and dependencies:
        # The Python code uses classes. In R, R6 classes are used to replicate class behavior.
        # The Python Bio module (Bio.Align, Bio.Align.substitution_matrices) is mapped to the Biostrings package in R for sequence alignment using the pairwiseAlignment function.
        # The BLOSUM62 substitution matrix is accessed via Biostrings::BLOSUM62.
        # The JSON functionality in Python (json.load) is provided by the jsonlite package in R (fromJSON).
# Type system differences:
        # Python’s dynamic typing is emulated with R’s flexible list structures.
        # The Python defaultdict is replaced by using lists and checking with is.null before initializing a new vector.
# Standard library equivalents:
        # os.path.join is replaced by file.path.
        # os.path.exists is replaced by file.exists.
        # The “with open(..., encoding)” idiom is replaced by readLines and concatenating the lines for JSON parsing.
# Error handling:
        # The original Python code does not implement explicit error handling, so no additional error handling was added in the R translation.
# Implementation decisions:
        # The swap of query_seq and ref_seq inside the get_alignment method is preserved exactly as per the original code.
        # For sequence alignment, default gap penalties (gapOpening = 11, gapExtension = 1) are used to mimic BLASTP as closely as possible.
        # Iterations and dictionary manipulations in Python are translated into for-loops and list manipulations in R.
        # All comments from the original Python code are preserved as comments in R.
# Potential issues and limitations:
        # The translation assumes that input sequences are provided as character strings.
        # The structure of lists (structure_list, fasta_list) in R is assumed to be similar to that in Python (i.e., lists of two-element lists/vectors).
        # The alignment parameters may not exactly replicate BLASTP’s behavior; adjustments might be needed for exact equivalence.
