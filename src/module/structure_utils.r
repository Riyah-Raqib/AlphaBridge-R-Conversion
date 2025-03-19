# ==============================
# Load Required Libraries
# ==============================
library(R6)
library(jsonlite)
library(ini)
library(reticulate)
library(Biostrings)
library(data.table)
library(stringr)
library(tools)
options(warn = -1)  # suppress warnings

# Import Python's numpy via reticulate for AF2 pickle file loading
np <- import("numpy", convert = TRUE)

# -------------------------
# Dummy stub implementations for external modules to preserve functionality
# -------------------------

# Stub for PDBPARSER with a get_ca_distances method
PDBPARSER <- R6Class("PDBPARSER",
    public = list(
        structure_path = NULL,
        initialize = function(structure_path) {
            self$structure_path <- structure_path
        },
        get_ca_distances = function() {
            # Dummy implementation: returns a 10x3 matrix of random numbers as coordinates
            coords <- matrix(runif(30), ncol = 3)
            # Compute pairwise Euclidean distances between rows
            return(as.matrix(dist(coords)))
        }
    )
)

# Stub for MMCIFPARSER with required methods
MMCIFPARSER <- R6Class("MMCIFPARSER",
    public = list(
        structure_path = NULL,
        initialize = function(structure_path) {
            self$structure_path <- structure_path
        },
        get_sequence_list = function() {
            # Dummy implementation: return a list with one record having id and seq
            return(list(list(id = "A", seq = "ACDEFGHIK")))
        },
        get_ca_distances = function() {
            # Dummy implementation: returns a 10x3 matrix of random numbers as coordinates
            coords <- matrix(runif(30), ncol = 3)
            return(as.matrix(dist(coords)))
        },
        get_polypeptide_chain_dict = function() {
            # Dummy implementation: returns a list with one chain "A"
            return(list(A = list(entity_type = "polypeptide(L)")))
        },
        get_coordinates = function() {
            # Dummy implementation: returns a nested list structure mimicking coordinates with plddt values
            return(list(A = list("1" = list(atom_id = list("N" = list(plddt = 90),
                                                                "CA" = list(plddt = 91))))))
        }
    )
)

# Stub for compare_protein_seq to mimic the alignment utility
compare_protein_seq <- function(structure_sequence_list, fasta_sequence_list) {
    extract_chain_dict <- function() {
        # Dummy implementation: map each fasta record id to itself if present in structure_sequence_list
        mapping <- list()
        for (rec in fasta_sequence_list) {
            mapping[[ rec[[1]] ]] <- rec[[1]]
        }
        return(mapping)
    }
    return(list(extract_chain_dict = extract_chain_dict))
}

# -------------------------
# End of stubs for external modules
# -------------------------

# Read configuration from config.ini using ini package
module_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[1]))
config <- read.ini(file.path(module_dir, "config.ini"))

WORKING_DIR <- config$DEFAULT$WORKING_DIR
AF3_DIR <- config$DEFAULT$AF3_DIR

# Define NumpyEncoder equivalent function to help with JSON encoding of matrices and arrays
NumpyEncoder <- function(x) {
    if (is.matrix(x) || is.array(x)) {
        return(as.list(x))
    }
    return(x)
}

FEATURE_MATRIX <- R6Class("FEATURE_MATRIX",
    public = list(
        job_id = NULL,
        alphafold_version = NULL,
        outdir = NULL,
        initialize = function(job_id, outdir = '', alphafold_version = 'AF2') {
            self$job_id <- job_id
            self$alphafold_version <- alphafold_version
            self$outdir <- outdir
        },
        get_sequence_info_path = function() {
            job_id <- self$job_id
            if (self$alphafold_version == 'AF2') {
                filepath <- sprintf("/DATA/zata/projects/interactions/data/request/fasta/%s.fasta", job_id)
            } else if (self$alphafold_version == 'AF3') {
                filepath <- self$get_feature_folder_path()
            }
            if (file.exists(filepath)) {
                return(filepath)
            } else {
                stop(sprintf("[%d] %s", errno(), strerror(errno(), filepath)))
            }
        },
        get_feature_folder_path = function() {
            job_id <- self$job_id
            if (self$alphafold_version == 'AF2') {
                filepath <- sprintf("/DATA/zata/projects/interactions/data/results/%s", job_id)
            } else if (self$alphafold_version == 'AF3') {
                filepath <- file.path(AF3_DIR, job_id)
            }
            if (file.exists(filepath)) {
                return(filepath)
            } else {
                stop(sprintf("[%d] %s", errno(), strerror(errno(), filepath)))
            }
        },
        extract_sequences = function(sequence_info_path) {
            if (self$alphafold_version == 'AF2') {
                # Using Biostrings to read fasta file
                recs <- readAAStringSet(sequence_info_path)
                rec_list <- list()
                for (i in seq_along(recs)) {
                    rec <- list(
                        id = names(recs)[i],
                        seq = as.character(recs[[i]])
                    )
                    rec_list[[length(rec_list) + 1]] <- rec
                }
            } else if (self$alphafold_version == 'AF3') {
                folder_path <- sequence_info_path
                job_request_pattern <- "job_request.*\\.json"
                files <- list.files(folder_path, pattern = job_request_pattern, full.names = TRUE)
                if (length(files) == 0) {
                    stop("No job_request json file found.")
                }
                job_request_all <- fromJSON(files[1])
                # In Python, job_request is set to the first element of the loaded json array
                job_request <- job_request_all[[1]]
                chains <- LETTERS
                chain_counter <- 1
                known_keys <- c("proteinChain", "ligand", "rnaSequence", "ion", "dnaSequence")
                seq_types_symbols <- list(proteinChain = "", rnaSequence = "RNA_", dnaSequence = "DNA_", ligand = "ligand_", ion = "ion_")
                rec_list <- list()
                for (macromolecule in job_request$sequences) {
                    key <- names(macromolecule)[1]
                    count <- macromolecule[[key]]$count
                    for (item in 1:count) {
                        if (key %in% c('proteinChain', 'dnaSequence', 'rnaSequence')) {
                            seq_obj <- macromolecule[[key]]$sequence
                            chain_id <- paste0(seq_types_symbols[[key]], chains[chain_counter])
                            chain_counter <- chain_counter + 1
                            rec <- list(id = chain_id, name = chain_id, description = chain_id, seq = seq_obj)
                            rec_list[[length(rec_list) + 1]] <- rec
                        }
                    }
                }
            }
            return(rec_list)
        },
        get_sequence_chain_tuple = function() {
            sequence_info_path <- self$get_sequence_info_path()
            sequence_list <- self$extract_sequences(sequence_info_path)
            sequence_chain_tuple <- lapply(sequence_list, function(rec) {
                return(list(rec$id, rec$seq))
            })
            return(sequence_chain_tuple)
        },
        reorder_sequences_by_token_list = function() {
            files <- self$find_feature_files()
            structure_path <- files[[1]]
            feature_path <- files[[2]]
            structure_sequence_list <- MMCIFPARSER$new(structure_path)$get_sequence_list()
            fasta_sequence_list <- self$get_sequence_chain_tuple()
            fasta_to_struct <- compare_protein_seq(structure_sequence_list, fasta_sequence_list)$extract_chain_dict()
            struct_to_fasta <- list()
            for (k in names(fasta_to_struct)) {
                struct_to_fasta[[ as.character(fasta_to_struct[[k]]) ]] <- k
            }
            feature_dict <- read_json_file(feature_path)
            token_chain_ids <- feature_dict$token_chain_ids
            unique_chains <- unique(token_chain_ids)
            unique_fasta_chains <- c()
            for (unique_chain in unique_chains) {
                if (unique_chain %in% names(struct_to_fasta)) {
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
            reorder_rec_list <- lapply(reorder_fasta_sequence_list, function(reorder_tuple) {
                return(list(id = reorder_tuple[[1]], name = reorder_tuple[[1]], description = reorder_tuple[[1]], seq = reorder_tuple[[2]]))
            })
            return(reorder_rec_list)
        },
        fasta_profiles = function() {
            list_fasta_name <- c()
            list_fasta_acclen <- c()
            list_fasta_len <- c()
            list_fasta_centerticks <- c()
            num_acc <- 0
            sequence_info_path <- self$get_sequence_info_path()
            if (self$alphafold_version == 'AF2') {
                fasta_sequences <- self$extract_sequences(sequence_info_path)
            } else if (self$alphafold_version == 'AF3') {
                fasta_sequences <- self$reorder_sequences_by_token_list()
            }
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
        find_feature_files = function() {
            if (self$alphafold_version == 'AF2') {
                feature_folder <- self$get_feature_folder_path()
                ranking_file <- file.path(feature_folder, 'ranking_debug.json')
                json_text <- readLines(ranking_file, warn = FALSE)
                ranking <- fromJSON(paste(json_text, collapse = "\n"))
                first_model <- ranking$order[1]
                feature_filename <- sprintf("result_%s.pkl", first_model)
                structure_path <- file.path(feature_folder, "ranked_0.pdb")
                feature_path <- file.path(feature_folder, feature_filename)
            } else if (self$alphafold_version == 'AF3') {
                folder_path <- self$get_feature_folder_path()
                json_files <- list.files(folder_path, pattern = "full_data_0\\.json", full.names = TRUE)
                cif_files <- list.files(folder_path, pattern = "model_0\\.cif", full.names = TRUE)
                if (length(json_files) == 0 || length(cif_files) == 0) {
                    stop("Required feature files not found in AF3 directory.")
                }
                feature_path <- json_files[1]
                structure_path <- cif_files[1]
            }
            return(list(structure_path, feature_path))
        },
        get_distance_matrix = function() {
            files <- self$find_feature_files()
            structure_path <- files[[1]]
            feature_path <- files[[2]]
            if (self$alphafold_version == 'AF2') {
                ca_distances <- PDBPARSER$new(structure_path)$get_ca_distances()
            } else if (self$alphafold_version == 'AF3') {
                ca_distances <- MMCIFPARSER$new(structure_path)$get_ca_distances()
            }
            distance_matrix <- as.matrix(dist(ca_distances))
            return(distance_matrix)
        },
        fix_matrix_size = function(structure_path, feature_dict) {
            pae <- as.matrix(feature_dict$pae)
            contact_probability <- as.matrix(feature_dict$contact_probs)
            token_chain_ids <- feature_dict$token_chain_ids
            unique_chains <- unique(token_chain_ids)
            chain_index_dict <- list()
            for (chain in unique_chains) {
                ii <- which(token_chain_ids == chain)
                chain_index_dict[[chain]] <- ii
            }
            polymer_chain_dict <- MMCIFPARSER$new(structure_path)$get_polypeptide_chain_dict()
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
            fix_size_pae <- pae[mask, mask, drop = FALSE]
            fix_size_contact_probability <- contact_probability[mask, mask, drop = FALSE]
            return(list(fix_size_pae, fix_size_contact_probability))
        },
        extract_plddt_per_residue = function(structure_path) {
            structure_coordinates <- MMCIFPARSER$new(structure_path)$get_coordinates()
            polymer_chain_dict <- MMCIFPARSER$new(structure_path)$get_polypeptide_chain_dict()
            data_list <- list()
            for (asym_id in names(structure_coordinates)) {
                for (seq_id in names(structure_coordinates[[asym_id]])) {
                    atom_ids <- names(structure_coordinates[[asym_id]][[seq_id]]$atom_id)
                    for (atom_id in atom_ids) {
                        if (polymer_chain_dict[[asym_id]]$entity_type %in% c('polypeptide(L)', 'polydeoxyribonucleotide', 'polyribonucleotide')) {
                            plddt <- as.numeric(structure_coordinates[[asym_id]][[seq_id]]$atom_id[[atom_id]]$plddt)
                            data_list[[length(data_list) + 1]] <- list(asym_id = asym_id, seq_id = seq_id, plddt = plddt)
                        }
                    }
                }
            }
            plddt_df <- rbindlist(data_list)
            plddt_per_residue_df <- aggregate(plddt ~ asym_id + seq_id, data = plddt_df, FUN = mean)
            residue_plddts <- plddt_per_residue_df$plddt
            return(residue_plddts)
        },
        get_feature_info = function() {
            files <- self$find_feature_files()
            structure_path <- files[[1]]
            feature_path <- files[[2]]
            if (self$alphafold_version == 'AF2') {
                distance_matrix <- self$get_distance_matrix()
                # Load pickle file using reticulate's numpy load
                feature_dict <- np$load(feature_path, allow_pickle = TRUE)
                # Convert to R list
                feature_dict <- py_to_r(feature_dict)
                pae <- feature_dict$predicted_aligned_error
                plddt <- feature_dict$plddt
                contact_probability <- list()
            } else if (self$alphafold_version == 'AF3') {
                folder_path <- self$get_feature_folder_path()
                feature_dict <- read_json_file(feature_path)
                plddt <- self$extract_plddt_per_residue(structure_path)
                distance_matrix <- self$get_distance_matrix()
                fix_size <- self$fix_matrix_size(structure_path, feature_dict)
                pae <- fix_size[[1]]
                contact_probability <- fix_size[[2]]
            }
            return(list(distance_matrix, pae, contact_probability, plddt))
        },
        extract_matrix_dict = function() {
            matrix_dict <- list()
            feat_info <- self$get_feature_info()
            distance_matrix <- feat_info[[1]]
            pae <- feat_info[[2]]
            contact_probability <- feat_info[[3]]
            plddt <- feat_info[[4]]
            symmetric_pae <- pae
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
                    if ((delta_index >= -2) && (delta_index <= 2)) {
                        plddt_matrix[i, j] <- 0
                    } else {
                        plddt_matrix[i, j] <- -1 * (((plddt[i] + plddt[j]) / 2) - 100)
                    }
                }
            }
            pae_plddt <- symmetric_pae + plddt_matrix / 3
            pae_plddt[pae_plddt > 32] <- 32
            if (self$alphafold_version == 'AF2') {
                log_distance <- ifelse(distance_matrix != 0, log10(distance_matrix), 0)
                log_pae <- ifelse(symmetric_pae != 0, log10(symmetric_pae), 0)
                confidance_matrix <- log_distance + log_pae
                contact_matrix <- distance_matrix
                modified_distance_matrix <- distance_matrix
                modified_distance_matrix[modified_distance_matrix > 40] <- 40
                mask_upper <- upper.tri(modified_distance_matrix, diag = TRUE)
                masked_contact_matrix <- modified_distance_matrix
                masked_contact_matrix[mask_upper] <- NA
                mask_lower <- lower.tri(symmetric_pae, diag = TRUE)
                masked_confidance_matrix <- pae
                masked_confidance_matrix[mask_lower] <- NA
            } else if (self$alphafold_version == 'AF3') {
                confidance_matrix <- pae_plddt
                contact_matrix <- contact_probability
                binary_contact <- contact_probability > 0.5
                mask_upper <- upper.tri(binary_contact, diag = TRUE)
                masked_contact_matrix <- binary_contact
                masked_contact_matrix[mask_upper] <- NA
                mask_lower <- lower.tri(pae_plddt, diag = TRUE)
                masked_confidance_matrix <- pae_plddt
                masked_confidance_matrix[mask_lower] <- NA
            }
            matrix_dict$pae <- pae
            matrix_dict$plddt <- plddt
            matrix_dict$plddt_matrix <- plddt_matrix
            matrix_dict$pae_plddt <- pae_plddt
            matrix_dict$symmetric_pae <- symmetric_pae
            matrix_dict$contact_matrix <- contact_matrix
            matrix_dict$confidance_matrix <- confidance_matrix
            matrix_dict$masked_confidance_matrix <- masked_confidance_matrix
            matrix_dict$masked_contact_matrix <- masked_contact_matrix
            return(matrix_dict)
        },
        print_matrix_dict = function() {
            matrix_dict <- self$extract_matrix_dict()
            matrix_dict$masked_confidance_matrix <- NULL
            matrix_dict$masked_contact_matrix <- NULL
            if (dir.exists(self$outdir)) {
                feature_object_path <- file.path(self$outdir, 'matrix_info.json')
                json_text <- toJSON(matrix_dict, pretty = TRUE, auto_unbox = TRUE, 
                                      force = TRUE, null = "null")
                write(json_text, file = feature_object_path)
            }
        }
    )
)

read_json_file <- function(json_file) {
    con <- file(json_file, "r")
    json_text <- readLines(con, warn = FALSE)
    close(con)
    file_data <- fromJSON(paste(json_text, collapse = "\n"))
    return(file_data)
}

# Helper functions to mimic errno and strerror in R for error messages
errno <- function() {
    return(2)
}

strerror <- function(err_no, filepath = "") {
    return(sprintf("No such file or directory: %s", filepath))
}

# Dependency Imports and Libraries:
# In this R translation, we load essential libraries: R6 (for class definitions), jsonlite (for JSON parsing), ini (for configuration parsing), reticulate (to load Python’s numpy for reading pickle files), Biostrings (for FASTA sequence handling), along with data.table, stringr, and tools for supporting functions.
# The warnings are suppressed using options(warn = -1) to mimic Python’s warnings.filterwarnings("ignore").
# External Module Stubs:
    # The Python code imports external modules (PDBPARSER, MMCIFPARSER, and compare_protein_seq). Since their full implementations are not provided, I created R6 stub classes for PDBPARSER and MMCIFPARSER with the required methods that return dummy data.
    # Similarly, compare_protein_seq is implemented as a simple function that returns an identity mapping.
    # These stubs ensure that every call in the original code is preserved and the translation is complete. However, in a production environment, these stubs must be replaced with proper implementations.
# Configuration Parsing:
    # The configuration file is read using the ini package. The module directory is determined using commandArgs and normalizePath. This may need adjustment depending on the deployment environment.
# Class and Method Translations:
    # All Python class methods are translated into methods of an R6 class named FEATURE_MATRIX while preserving method names and functionality.
    # Python’s string formatting (f-strings) is replaced with sprintf.
    # File path operations use file.path and related R functions.
    # Error handling in Python via exceptions is replaced with the R stop() function.
# Sequence Handling:
    # For extracting sequences, Biostrings is used to read FASTA files when alphafold_version is "AF2".
    # When alphafold_version is "AF3", the code uses list.files to perform glob-style matching, and the sequence records are manually created as lists to mimic Bio.SeqRecord. A counter based on LETTERS is used to simulate chain identification.
# Matrix and Array Operations:
    #Numpy operations are translated to R’s matrix operations. For pairwise distances, the as.matrix(dist(...)) function is used.
    #The symmetric pae matrix and plddt matrix operations are implemented using nested for loops to closely mimic the Python code.
    #Masked arrays are simulated by setting masked elements to NA.
#JSON and Pickle Handling:
    #JSON reading and writing are done using jsonlite.
    #For AF2, Python’s pickle file is loaded using reticulate to call numpy’s load function. The py_to_r function converts the loaded data to an R list.
    #The custom NumpyEncoder class in Python is replaced with the NumpyEncoder function for converting matrices/arrays during JSON conversion, although jsonlite handles most conversions automatically.
    #Preservation of Code Structure:
    #Every line and comment from the original Python code has been translated and preserved as closely as possible.
    #All variables, methods, and class names are maintained, ensuring an exact translation.
#Potential Issues and Limitations:
    #The stub implementations for external modules (PDBPARSER, MMCIFPARSER, compare_protein_seq) provide dummy data. These will need proper implementations for real data processing.
    #The method of obtaining the module directory in R may require adjustments depending on how the script is executed.
    #The handling of pickle files via reticulate assumes that a working Python environment is available.
    #The simulation of masked arrays (using NA) is simplistic compared to NumPy’s masked arrays.
