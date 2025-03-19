# ==============================
# Load Required Libraries
# ==============================
library(R6)
library(jsonlite)
library(tools)    # for file_path_sans_ext if needed

# Suppress warnings, equivalent to Python warnings.filterwarnings("ignore")
options(warn = -1)

# Dummy definitions for external dependencies
# These are provided to preserve the structure of the original code.
# Users must implement or source the actual implementations.
PDBPARSER <- function(filepath) {
  # Dummy implementation for PDBPARSER
  DummyPDB <- R6Class("DummyPDB",
    public = list(
      get_sequence_list = function() { },
      get_polypeptide_chain_dict = function() { },
      get_coordinates = function() { },
      get_ca_distances = function() { }
    )
  )
  DummyPDB$new()
}

MMCIFPARSER <- function(filepath) {
  # Dummy implementation for MMCIFPARSER
  DummyMMCIF <- R6Class("DummyMMCIF",
    public = list(
      get_sequence_list = function() { },
      get_polypeptide_chain_dict = function() { },
      get_coordinates = function() { },
      get_ca_distances = function() { }
    )
  )
  DummyMMCIF$new()
}

compare_protein_seq <- function(structure_sequence_list, fasta_sequence_list) {
  # Dummy implementation for compare_protein_seq
  DummyCompare <- R6Class("DummyCompare",
    public = list(
      extract_chain_dict = function() { list() }
    )
  )
  DummyCompare$new()
}

# Define module_dir equivalent
module_dir <- dirname(normalizePath(getwd()))

# NumpyEncoder class equivalent in R.
NumpyEncoder <- R6Class("NumpyEncoder",
  public = list(
    default = function(obj) {
      # if obj is a numeric vector or matrix, convert to list
      if (is.matrix(obj) || is.numeric(obj)) {
        return(as.list(obj))
      }
    }
  )
)

FEATURE_MATRIX <- R6Class("FEATURE_MATRIX",
  public = list(
    in_dir = NULL,
    
    initialize = function(in_dir) {
      self$in_dir <- in_dir
    },
    
    check_if_path_exist = function(filepath) {
      # Check if the file or directory exists
      if (file.exists(filepath)) {
        return(TRUE)
      } else {
        stop(paste("FileNotFoundError:", filepath, "does not exist"))
      }
    },
    
    fasta_profiles = function(fasta_sequences) {
      list_fasta_name <- c()
      list_fasta_acclen <- c()
      list_fasta_len <- c()
      list_fasta_files <- c()
      list_fasta_centerticks <- c()
      num_acc <- 0
      
      for (i in seq_along(fasta_sequences)) {
        fasta <- fasta_sequences[[i]]
        name <- fasta$id
        sequence <- as.character(fasta$seq)
        list_fasta_name <- c(list_fasta_name, name)
        num_acc <- num_acc + nchar(sequence)
        if (length(list_fasta_acclen) == 0) {
          center_tick <- as.integer((num_acc - 0) / 2)
        } else {
          center_tick <- as.integer((num_acc - tail(list_fasta_acclen, n=1)) / 2 + tail(list_fasta_acclen, n=1))
        }
        list_fasta_centerticks <- c(list_fasta_centerticks, center_tick)
        list_fasta_acclen <- c(list_fasta_acclen, num_acc)
        list_fasta_len <- c(list_fasta_len, nchar(sequence))
      }
      
      return(list(list_fasta_name, list_fasta_acclen, list_fasta_centerticks, list_fasta_len))
    },
    
    get_sequence_chain_tuple = function(sequence_list) {
      sequence_chain_tuple <- lapply(sequence_list, function(rec) {
        list(rec$id, rec$seq)
      })
      return(sequence_chain_tuple)
    },
    
    get_distance_matrix = function(ca_distances) {
      # Compute the pairwise Euclidean distances.
      # Assumes ca_distances is a numeric matrix/data frame with observations in rows.
      distance_matrix <- as.matrix(dist(ca_distances))
      return(distance_matrix)
    },
    
    get_pae_plddt_matrix = function(pae, plddt) {
      symmetric_pae <- pae  # copy
      
      # Adjust symmetric_pae values to be symmetric by taking the maximum in the two directions
      for (i in 1:nrow(symmetric_pae)) {
        for (j in 1:ncol(symmetric_pae)) {
          if (symmetric_pae[i, j] < symmetric_pae[j, i]) {
            symmetric_pae[i, j] <- symmetric_pae[j, i]
          }
        }
      }
      
      plddt_matrix <- matrix(0, nrow = length(plddt), ncol = length(plddt))
      
      for (i in 1:nrow(plddt_matrix)) {
        for (j in 1:ncol(plddt_matrix)) {
          delta_index <- i - j
          if (delta_index >= -2 && delta_index <= 2) {
            plddt_matrix[i, j] <- 0
          } else {
            plddt_matrix[i, j] <- -1 * (((plddt[i] + plddt[j]) / 2) - 100)
          }
        }
      }
      
      pae_plddt <- symmetric_pae + plddt_matrix / 3
      
      # Cap values in pae_plddt at 32
      pae_plddt[pae_plddt > 32] <- 32
      
      return(list(symmetric_pae, pae_plddt, plddt_matrix))
    },
    
    get_feature_matrix_dict = function(pae, plddt, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidance_matrix, masked_confidance_matrix, masked_contact_matrix) {
      matrix_dict <- list()
      matrix_dict[["pae"]] <- pae
      matrix_dict[["plddt"]] <- plddt
      matrix_dict[["plddt_matrix"]] <- plddt_matrix
      matrix_dict[["pae_plddt"]] <- pae_plddt
      matrix_dict[["symmetric_pae"]] <- symmetric_pae
      matrix_dict[["contact_matrix"]] <- contact_matrix
      matrix_dict[["confidence_matrix"]] <- confidance_matrix
      matrix_dict[["masked_confidence_matrix"]] <- masked_confidance_matrix
      matrix_dict[["masked_contact_matrix"]] <- masked_contact_matrix
      
      return(matrix_dict)
    },
    
    get_scores_dict = function(scores_list, list_fasta_files) {
      scores_dict <- list()
      initial_acc <- 1  # R indexing starts at 1
      
      temp <- list_fasta_files
      list_fasta_name <- temp[[1]]
      list_fasta_acclen <- temp[[2]]
      
      for (k in seq_along(list_fasta_name)) {
        fasta_name <- list_fasta_name[k]
        fasta_acclen <- list_fasta_acclen[k]
        if (is.null(scores_dict[[fasta_name]])) {
          scores_dict[[fasta_name]] <- scores_list[initial_acc:fasta_acclen]
        }
        initial_acc <- fasta_acclen + 1
      }
      
      return(scores_dict)
    },
    
    print_matrix_dict = function(matrix_dict) {
      excluded_keys <- c('pae','plddt','symmetric_pae','pae_plddt','masked_confidance_matrix', 'masked_contact_matrix')
      tmp_matrix <- matrix_dict[setdiff(names(matrix_dict), excluded_keys)]
      
      if (file.exists(self$in_dir)) {
        feature_object_path <- file.path(self$in_dir, 'matrix_info.json')
        # Write JSON to file using the NumpyEncoder equivalent
        encoder <- NumpyEncoder$new()
        json_text <- toJSON(tmp_matrix, auto_unbox = TRUE)
        write(json_text, file = feature_object_path)
      }
    }
  )
)

CCM_AF3 <- R6Class("CCM_AF3",
  inherit = FEATURE_MATRIX,
  public = list(
    initialize = function(in_dir) {
      super$initialize(in_dir)
    },
    
    extract_feature_filepath = function() {
      folder_path <- self$in_dir
      
      if (self$check_if_path_exist(folder_path)) {
        feature_files <- list.files(folder_path, pattern = "full_data_0\\.json$", full.names = TRUE)
        structure_files <- list.files(folder_path, pattern = "model_0\\.cif$", full.names = TRUE)
        job_request_files <- list.files(folder_path, pattern = "job_request.*\\.json$", full.names = TRUE)
        
        feature_path <- feature_files[1]
        structure_path <- structure_files[1]
        job_request_path <- job_request_files[1]
        
        return(list(feature_path, structure_path, job_request_path))
      }
    },
    
    extract_sequences = function(job_request_path) {
      request_file <- fromJSON(job_request_path)[[1]]
      
      chains <- LETTERS  # Generator for A-Z
      known_keys <- c("proteinChain", "ligand", "rnaSequence", "ion", "dnaSequence")
      seq_types_symbols <- list("proteinChain" = "", "rnaSequence" = "RNA_", "dnaSequence" = "DNA_", "ligand" = "ligand_", "ion" = "ion_")
      
      rec_list <- list()
      chain_index <- 1
      
      for (macromolecule in request_file$sequences) {
        key <- names(macromolecule)[1]
        count <- macromolecule[[key]]$count
        if (key %in% c('proteinChain', 'dnaSequence', 'rnaSequence')) {
          for (item in 1:count) {
            seq <- macromolecule[[key]]$sequence
            chain_id <- paste0(seq_types_symbols[[key]], chains[chain_index])
            chain_index <- chain_index + 1
            rec <- list(id = chain_id, seq = seq, name = chain_id, description = chain_id)
            rec_list[[length(rec_list) + 1]] <- rec
          }
        }
      }
      
      return(rec_list)
    },
    
    extract_sequence_info = function() {
      paths <- self$extract_feature_filepath()
      feature_path <- paths[[1]]
      structure_path <- paths[[2]]
      job_request_path <- paths[[3]]
      
      structure <- MMCIFPARSER(structure_path)
      feature_dict <- read_json_file(feature_path)
      
      structure_sequence_list <- structure$get_sequence_list()
      polymer_chain_dict <- structure$get_polypeptide_chain_dict()
      
      tmp_rec_list <- self$extract_sequences(job_request_path)
      rec_list <- self$reorder_sequences_by_token_list(tmp_rec_list, structure_sequence_list, feature_dict)
      chain_tuple <- self$get_sequence_chain_tuple(rec_list)
      list_sequence_info <- self$fasta_profiles(rec_list)
      
      return(list(list_sequence_info, chain_tuple, structure_sequence_list, polymer_chain_dict))
    },
    
    reorder_sequences_by_token_list = function(tmp_rec_list, structure_sequence_list, feature_dict) {
      fasta_sequence_list <- self$get_sequence_chain_tuple(tmp_rec_list)
      
      comp_obj <- compare_protein_seq(structure_sequence_list, fasta_sequence_list)
      fasta_to_struct <- comp_obj$extract_chain_dict()
      
      # Invert the dictionary: keys become values and vice versa
      struct_to_fasta <- list()
      for (k in names(fasta_to_struct)) {
        struct_to_fasta[[ fasta_to_struct[[k]] ]] <- k
      }
      
      token_chain_ids <- feature_dict$token_chain_ids
      unique_chains <- unique(token_chain_ids)
      
      unique_fasta_chains <- c()
      for (unique_chain in unique_chains) {
        if (!is.null(struct_to_fasta[[unique_chain]])) {
          unique_fasta_chains <- c(unique_fasta_chains, struct_to_fasta[[unique_chain]])
        }
      }
      
      reorder_fasta_sequence_list <- list()
      for (x in unique_fasta_chains) {
        for (tuple in fasta_sequence_list) {
          if (tuple[[1]] == x) {
            reorder_fasta_sequence_list[[length(reorder_fasta_sequence_list) + 1]] <- tuple
          }
        }
      }
      
      reorder_rec_list <- list()
      for (reorder_tuple in reorder_fasta_sequence_list) {
        rec <- list(id = reorder_tuple[[1]], seq = reorder_tuple[[2]], name = reorder_tuple[[1]], description = reorder_tuple[[1]])
        reorder_rec_list[[length(reorder_rec_list) + 1]] <- rec
      }
      
      return(reorder_rec_list)
    },
    
    extract_plddt_per_residue = function(structure) {
      structure_coordinates <- structure$get_coordinates()
      polymer_chain_dict <- structure$get_polypeptide_chain_dict()
      
      data <- list()
      
      for (asym_id in names(structure_coordinates)) {
        for (seq_id in names(structure_coordinates[[asym_id]])) {
          atom_ids <- names(structure_coordinates[[asym_id]][[seq_id]]$atom_id)
          for (atom_id in atom_ids) {
            if (polymer_chain_dict[[asym_id]]$entity_type %in% c('polypeptide(L)', 'polydeoxyribonucleotide', 'polyribonucleotide')) {
              plddt <- as.numeric(structure_coordinates[[asym_id]][[seq_id]]$atom_id[[atom_id]]$plddt)
              data[[length(data) + 1]] <- list(asym_id = asym_id, seq_id = seq_id, plddt = plddt)
            }
          }
        }
      }
      
      # Convert list of records to a data frame
      if (length(data) > 0) {
        plddt_df <- do.call(rbind, lapply(data, as.data.frame))
        # Group by asym_id and seq_id and compute mean of plddt
        plddt_per_residue_df <- aggregate(plddt ~ asym_id + seq_id, data = plddt_df, FUN = mean)
        residue_plddts <- as.numeric(plddt_per_residue_df$plddt)
      } else {
        residue_plddts <- numeric(0)
      }
      
      return(residue_plddts)
    },
    
    fix_matrix_size = function(structure, feature_dict) {
      pae <- as.matrix(feature_dict$pae)
      contact_probability <- as.matrix(feature_dict$contact_probs)
      
      token_chain_ids <- feature_dict$token_chain_ids
      unique_chains <- unique(token_chain_ids)
      chain_index_dict <- list()
      
      for (chain in unique_chains) {
        ii <- which(token_chain_ids == chain)
        chain_index_dict[[chain]] <- ii
      }
      
      # get chains that are a polypeptide
      polymer_chain_dict <- structure$get_polypeptide_chain_dict()
      polypeptide_chain_list <- c()
      for (chain in names(polymer_chain_dict)) {
        if (polymer_chain_dict[[chain]]$entity_type %in% c('polypeptide(L)', 'polydeoxyribonucleotide', 'polyribonucleotide')) {
          polypeptide_chain_list <- c(polypeptide_chain_list, chain)
        }
      }
      
      non_polypeptide_chain_list <- setdiff(unique_chains, polypeptide_chain_list)
      
      remove_idx_list <- c()
      
      for (non_polypeptide_chain in non_polypeptide_chain_list) {
        remove_idx_list <- c(remove_idx_list, chain_index_dict[[non_polypeptide_chain]])
      }
      
      mask <- rep(TRUE, nrow(pae))
      mask[remove_idx_list] <- FALSE
      
      fix_size_pae <- pae[mask, mask, drop=FALSE]
      fix_size_contact_probability <- contact_probability[mask, mask, drop=FALSE]
      
      return(list(fix_size_pae, fix_size_contact_probability))
    },
    
    get_feature_info = function() {
      paths <- self$extract_feature_filepath()
      feature_path <- paths[[1]]
      structure_path <- paths[[2]]
      job_request_path <- paths[[3]]
      
      structure <- MMCIFPARSER(structure_path)
      
      feature_dict <- read_json_file(feature_path)
      
      plddt <- self$extract_plddt_per_residue(structure)
      
      distance_matrix <- self$get_distance_matrix(structure$get_ca_distances())
      
      fixed_matrices <- self$fix_matrix_size(structure, feature_dict)
      pae <- fixed_matrices[[1]]
      contact_probability <- fixed_matrices[[2]]
      
      return(list(distance_matrix, pae, contact_probability, plddt))
    },
    
    extract_matrix_dict = function() {
      feat_info <- self$get_feature_info()
      distance_matrix <- feat_info[[1]]
      pae <- feat_info[[2]]
      contact_probability <- feat_info[[3]]
      plddt <- feat_info[[4]]
      
      pae_out <- self$get_pae_plddt_matrix(pae, plddt)
      symmetric_pae <- pae_out[[1]]
      pae_plddt <- pae_out[[2]]
      plddt_matrix <- pae_out[[3]]
      
      confidance_matrix <- pae_plddt
      contact_matrix <- contact_probability
      binary_contact <- contact_matrix > 0.5
      
      # Create mask for upper triangle including diagonal
      mask_upper <- upper.tri(binary_contact, diag = TRUE)
      masked_contact_matrix <- binary_contact
      masked_contact_matrix[mask_upper] <- NA  # Set masked entries to NA
      
      # Create mask for lower triangle including diagonal
      mask_lower <- lower.tri(pae_plddt, diag = TRUE)
      masked_confidance_matrix <- pae_plddt
      masked_confidance_matrix[mask_lower] <- NA  # Set masked entries to NA
      
      matrix_dict <- self$get_feature_matrix_dict(pae, plddt, plddt_matrix, pae_plddt, symmetric_pae, contact_matrix, confidance_matrix, masked_confidance_matrix, masked_contact_matrix)
      
      return(matrix_dict)
    }
  )
)

read_json_file <- function(json_file) {
  con <- file(json_file, "r")
  file_content <- readLines(con, warn = FALSE)
  close(con)
  file_parsed <- fromJSON(paste(file_content, collapse = ""))
  return(file_parsed)
}

# Language-specific feature adaptations:
    # The Python classes have been translated to R6 classes in R to closely mimic the object-oriented structure.
    # All Python methods are implemented as R6 public methods with the same names.
    # The nested loops and matrix manipulations have been translated using R’s for-loops and matrix indexing (which is 1-indexed).
    # JSON operations are handled using the jsonlite package.
    # File path handling uses base R functions such as file.exists() and list.files() to mimic the functionality of pathlib.
# Type System Differences:
    # Python’s dynamic typing is mimicked with R’s flexible types. Where Python uses lists for tuples, R uses lists.
    # Numpy arrays and matrix operations have been mapped to R matrices.
# Standard Library Equivalents:
    # The Python json module is replaced with jsonlite.
    # The sklearn.metrics.pairwise.pairwise_distances is replaced by R’s as.matrix(dist(...)) function.
# Error Handling Approaches:
    # Python’s FileNotFoundError is replaced with stop() in R.
# External Dependencies and Dummy Implementations:
    # Dummy implementations for PDBPARSER, MMCIFPARSER, and compare_protein_seq have been provided since their implementations are assumed external.
    # These stubs must be replaced with actual implementations for real use.
# Masked Array Functionality:
    # Python’s np.ma.array has no direct equivalent in base R; therefore, the masked arrays have been simulated by setting masked values to NA.
# Comments and Formatting:
    # All comments from the original Python code have been preserved and translated to R-style comments using the # symbol.
    # The overall formatting, indentation, and structure of the code have been maintained to reflect the original layout.
# Potential Limitations:
    # The dummy implementations for external modules (PDBPARSER, MMCIFPARSER, etc.) are placeholders.
    # Some functions (especially those interacting with biological data structures) assume that the corresponding methods (e.g., get_sequence_list, get_coordinates) are defined in the external modules.
    # The translation assumes that the original logic and data structures in Python will map directly to R data structures without modification.
