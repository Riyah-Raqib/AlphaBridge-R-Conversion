# ==============================
# Load Required Libraries
# ==============================
library(R6)         # Modern object-oriented programming in R
library(Biostrings) # Handling protein and nucleotide sequences
library(stringr)    # String manipulation functions
library(dplyr)      # Data manipulation (e.g., filtering, summarisation)
library(tidyr)      # Data transformation (e.g., reshaping)
library(reticulate) # Interfacing with Python if needed

# ==============================
# Define Constants
# ==============================
# Mapping of 3-letter <-> 1-letter amino acid codes
protein_letters_3to1_extended <- c(
  "ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D", "CYS" = "C",
  "GLN" = "Q", "GLU" = "E", "GLY" = "G", "HIS" = "H", "ILE" = "I",
  "LEU" = "L", "LYS" = "K", "MET" = "M", "PHE" = "F", "PRO" = "P",
  "SER" = "S", "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V"
)

# Convert 3-letter codes to uppercase for consistency
upper_protein_letters_3to1 <- lapply(protein_letters_3to1_extended, toupper)

# ==============================
# General File Handling Function
# ==============================
# Read a file with checks for existence and supported file types
# Supports .gz (compressed), .pdb, .cif, and .hssp files
read_file_with_checks <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File not found - %s", filepath))
  }

  supported_extensions <- c("gz", "pdb", "cif", "hssp")
  ext <- tools::file_ext(filepath)

  if(!ext %in% supported_extensions) stop("Error: Unsupported file type")

  reader <- if (ext == "gz") gzfile else file()
  return(readLines(reader(filepath)))
}

# ==============================
# R6 Classes
# ==============================
# Msa Class -> Multiple sequence alignment
Msa <- R6Class("Msa",
  public = list(
    sequences = NULL,   #Charactor vector of sequences
    deletion_matrix = NULL,   #List of deletion matrices
    descriptions = NULL,    #Character vector of sequence descriptions

    # Constructor -> Initialise the MSA object
    initialise = function(sequences, deletion_matrix, descriptions) {
      if (!(length(sequences) == length(deletion_matrix) && length(sequences) == length(descriptions))) {
        stop("All fields for an MSA must have the same length.")
      }
      self$sequences <- sequences
      self$deletion_matrix <- deletion_matrix
      self$descriptions <- descriptions
    },

    # Truncate MSA to a specified number of sequences
    truncate = function(max_seqs) {
      return(Msa$new(
        sequences = self$sequences[1:max_seqs],
        deletion_matrix = self$deletion_matrix[1:max_seqs],
        descriptions = self$descriptions[1:max_seqs]
      ))
    }
  )
)

# TemplateHit Class -> Template hit in a sequence alignment
TemplateHit <- R6Class("TemplateHit",
  public = list(
    index = NULL,   # Integer index of the template hit
    name = NULL,    # Name of the template
    aligned_cols = NULL,    # Number of aligned columns
    sum_probs = NULL,   # Sum of probabilities for the alignment
    query = NULL,   # Query sentence
    hit_sequence = NULL,    # Hit sequence
    indices_query = NULL,   # Indices of the query sequence
    indices_hit = NULL,   # Indices of the hit sequence

    # Constructor -> Initialise the TemplateHit object
    initialise <- function(index, name, aligned_cols, sum_probs, query, hit_sequence, indices_query, indices_hit) {
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

# MMCIFPARSER Class -> Parse .cif files to extract structural information
MMCIFPARSER <- R6Class("MMCIFPARSER",
  public = list(
    filepath = NULL,    # Path to .cif file

    # Constructor -> Initialise the parser with a file path
    initialise <- function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

    # Read and return structure from .cif file
    get_structure_cif = function() {
      structure <- readMMCIF(self$filepath)
      return(structure)
    },

    # Extract atomic coordinates from .cif file
    get_coordinates = function() {
      structure <- readMMCIF(self$filepath)
      atoms <- structure$atom_site[structure$atom_site$group_PDB == "ATOM", ]
      coordinates_dict <- lapply(unique(atoms$label_asym_id), function(chain) {
        atoms[atoms$label_asym_id == chain, ]
      })
      names(coordinates_dict) <- unique(atoms$label_asym_id)
      return(coordinates_dict)
    }
  )
)

# PDBPARSER  -> Read and process /pdb files
PDBPARSER <- R6Class("PDBPARSER",
  public = list(
    filepath = NULL,    # Path to .pdb file

    # Constructor -> Initialise the parser with a file path
    initialise <- function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

      # Get atom coordinates from .pdb file
      get_coordinates = function() {
      lines <- self$get_lines_from_pdb_file()
      atom_lines <- grep("^ATOM", self$get_lines_from_pdb_file, value = TRUE)
      coordinates <- do.call(rbind, strsplit(substr(atom_lines, 31, 54), " "))
      coordinates <- matrix(as.numeric(coords), ncol = 3, byrow = TRUE)
      }
      base::return(coordinates)
    }

    # Read and return lines from .pdb file
    get_lines_from_pdb_file = function() {
      return(read_file_with_checks(self$filepath))
    }
  )
)

# HSSPPARSER Class -> Extract conservation information from .hssp files
HSSPPARSER <- R6Class("HSSPPARSER",
  public = list(
    filepath = NULL,    # Path to .hssp file

    # Constructor -> Initialise the parser with a file path
    initialise <- function(filepath) {
      if (!file.exists(filepath)) {
        stop(sprintf("File not found: %s", filepath))
      }
      self$filepath <- filepath
    },

    # Extract Shannon entropy, relative Shannon entropy, and sequence variability
    extract_conservation_info = function() {
      hssp_data <- readLines(self$filepath)

      matched_lines <- grep("^#=GF (PR|RI)", hssp_data, value = TRUE)

      # Extract relevant values
      entropy <- as.numeric(substr(matched_lines[grepl("PR", matched_lines)], 121, 125))
      relative_entropy <- entropy / log(20) * 100
      variability <- as.numeric(substr(matched_lines[grepl("RI", matched_lines)], 55, 57))

      return(list(
        shannon_entropy = shannon_entropy,
        relative_shannon_entropy = relative_shannon_entropy,
        sequence_variability = sequence_variability
      ))
    }
  )
)
