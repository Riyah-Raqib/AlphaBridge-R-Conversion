# ==============================
# Load Required Libraries
# ==============================
library(R6)         # Modern OOP implementation in R
library(Biostrings) # Handling protein and nucleotide sequences
library(stringr)    # String manipulation functions
library(dplyr)      # Data manipulation (e.g., filtering, summarisation)
library(tidyr)      # Data transformation (e.g., reshaping)
library(reticulate) # Interfacing with Python if needed

# ==============================
# Define Constants
# ==============================
protein_letters_3to1_extended <- list(
  "ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D", "CYS" = "C",
  "GLN" = "Q", "GLU" = "E", "GLY" = "G", "HIS" = "H", "ILE" = "I",
  "LEU" = "L", "LYS" = "K", "MET" = "M", "PHE" = "F", "PRO" = "P",
  "SER" = "S", "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V"
)

upper_protein_letters_3to1 <- lapply(protein_letters_3to1_extended, toupper)

# ==============================
# General File Handling Function
# ==============================
read_file_with_checks <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File not found - %s", filepath))
  }
  
  ext <- tools::file_ext(filepath)
  if (ext == "gz") {
    return(readLines(gzfile(filepath)))
  } else if (ext %in% c("pdb", "cif", "hssp")) {
    return(readLines(filepath))
  } else {
    stop("Error: Unsupported file type")
  }
}

# ==============================
# R6 Classes
# ==============================
# Msa Class
Msa <- R6Class("Msa",
  public = list(
    sequences = NULL,
    deletion_matrix = NULL,
    descriptions = NULL,
    
    initialize = function(sequences, deletion_matrix, descriptions) {
      if (!(length(sequences) == length(deletion_matrix) && length(sequences) == length(descriptions))) {
        stop("All fields for an MSA must have the same length.")
      }
      self$sequences <- sequences
      self$deletion_matrix <- deletion_matrix
      self$descriptions <- descriptions
    },
    
    truncate = function(max_seqs) {
      return(Msa$new(
        sequences = self$sequences[1:max_seqs],
        deletion_matrix = self$deletion_matrix[1:max_seqs],
        descriptions = self$descriptions[1:max_seqs]
      ))
    }
  )
)

# TemplateHit Class
TemplateHit <- R6Class("TemplateHit",
  public = list(
    index = NULL,
    name = NULL,
    aligned_cols = NULL,
    sum_probs = NULL,
    query = NULL,
    hit_sequence = NULL,
    indices_query = NULL,
    indices_hit = NULL,
    
    initialize = function(index, name, aligned_cols, sum_probs, query, hit_sequence, indices_query, indices_hit) {
      self$index <- index
      self$name <- name
      self$aligned_cols <- aligned_cols
      self$sum_probs <- sum_probs
      self$query <- query
      self$hit_sequence <- hit_sequence
      self$indices_query <- indices_query
      self$indices_hit <- indices_hit
    }
  )
)

# MMCIFPARSER Class
MMCIFPARSER <- R6Class("MMCIFPARSER",
  public = list(
    filepath = NULL,

    initialize = function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

    get_structure_cif = function() {
      structure <- readMMCIF(self$filepath)
      return(structure)
    },

    get_coordinates = function() {
      structure <- readMMCIF(self$filepath)
      atoms <- structure$atom_site[structure$atom_site$group_PDB == "ATOM", ]
      coordinates_dict <- split(atoms, atoms$label_asym_id)
      return(coordinates_dict)
    }
  )
)

# PDBPARSER Class
PDBPARSER <- R6Class("PDBPARSER",
  public = list(
    filepath = NULL,

    initialize = function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

    get_lines_from_pdb_file = function() {
      return(read_file_with_checks(self$filepath))
    }
  )
)

# HSSPPARSER Class
HSSPPARSER <- R6Class("HSSPPARSER",
  public = list(
    filepath = NULL,

    initialize = function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

    extract_conservation_info = function() {
      hssp_data <- readLines(self$filepath)
      
      shannon_entropy <- numeric()
      relative_shannon_entropy <- numeric()
      sequence_variability <- numeric()
      
      for (line in hssp_data) {
        if (startsWith(line, "#=GF PR")) {
          entropy <- as.numeric(substr(line, 121, 125))
          relative_entropy <- entropy / log(20) * 100
          shannon_entropy <- c(shannon_entropy, entropy)
          relative_shannon_entropy <- c(relative_shannon_entropy, relative_entropy)
        } else if (startsWith(line, "#=GF RI")) {
          sequence_variability <- c(sequence_variability, as.numeric(substr(line, 55, 57)))
        }
      }
      
      return(list(
        shannon_entropy = shannon_entropy,
        relative_shannon_entropy = relative_shannon_entropy,
        sequence_variability = sequence_variability
      ))
    }
  )
)
