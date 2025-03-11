# Required dependencies and imports
library(optparse)
library(tools)

# Source external modules (translated from Python "from src.module..." imports)
source("src/module/confidance_contact_matrix.R")    # Provides CCM_AF3
source("src/module/alingment_utils.R")               # Provides compare_protein_seq
source("src/module/domain_clustering.R")             # Provides domain_clustering
source("src/module/parsers.R")                       # Provides MMCIFPARSER, HSSPPARSER, alphafold_msa
# source("src/module/conservation_score.R")         # CONSERVATION_SCORE (commented out in original)
source("src/module/interface_identification.R")      # Provides interface_identification
source("src/module/ribbon_diagram.R")                # Provides RIBBON_DIAGRAM

# working_dir equivalent in R. In Python: obtain current file's directory.
# In R, we use getwd() as an approximation since __file__ is not available.
working_dir <- getwd()

parse_args <- function() {
    #####################
    # START CODING HERE #
    #####################
    # Implement a simple argument parser (WITH help documentation!) that parses
    # the information needed by main() from commandline.
    
    option_list <- list(
        make_option(c("-i"), dest = "in_dir", type = "character", default = "",
                    help = "path to a directory where input folder are stored"),
        make_option(c("-c"), dest = "config_dir", type = "character", default = "",
                    help = "path to a directory where input folder are stored"),
        make_option(c("-m", "--mode"), dest = "mode", type = "character", default = "AF3",
                    help = "output from different AlphaFold Version. Options: AF3, AF2, ColabFold"),
        make_option(c("-t", "--threshold"), dest = "contact_threshold", type = "double", default = 0.7,
                    help = "contact threshold to detect a contact-link in the contact_proability matrix")
    )
    
    parser <- OptionParser(option_list = option_list)
    
    # Check for no arguments provided and print help if so
    args_cmd <- commandArgs(trailingOnly = TRUE)
    if (length(args_cmd) == 0) {
        print_help(parser)
        stop()
    }
    
    args <- parse_args(parser)
    
    # Apply restricted_float conversion on contact_threshold
    args$contact_threshold <- restricted_float(args$contact_threshold)
    
    return(args)
}

restricted_float <- function(x) {
    # Try to convert x to numeric
    x_numeric <- as.numeric(x)
    if (is.na(x_numeric)) {
        stop(sprintf("'%s' not a floating-point literal", x))
    }
    
    if (x_numeric < 0.0 || x_numeric > 1.0) {
        stop(sprintf("'%s' not in range [0.0, 1.0]", x))
    }
    return(x_numeric)
}

write_dataframe <- function(df, filename, outdir_path) {
    filepath <- file.path(outdir_path, filename)
    write.csv(df, sprintf("%s.csv", filepath), row.names = FALSE)
}

define_interfaces <- function(in_dir, mode, contact_threshold) {
    outdir <- file.path(in_dir, "AlphaBridge")
    
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
    
    if (mode == "AF3") {
        FEATURE_OBJECT <- CCM_AF3(in_dir)
        
        temp <- FEATURE_OBJECT$extract_feature_filepath()
        feature_path <- temp$feature_path
        structure_path <- temp$structure_path
        job_request_path <- temp$job_request_path
        
        temp_seq <- FEATURE_OBJECT$extract_sequence_info()
        list_sequence_info <- temp_seq$list_sequence_info
        rec_sequence_list <- temp_seq$rec_sequence_list
        structure_sequence_list <- temp_seq$structure_sequence_list
        polymer_chain_dict <- temp_seq$polymer_chain_dict
    } else {
        stop("Output from AF2 or ColabFold not implemented yet")
    }
    
    chain_dict <- compare_protein_seq(structure_sequence_list, rec_sequence_list)$extract_chain_dict()
    
    matrix_dict <- FEATURE_OBJECT$extract_matrix_dict()
    
    plddt_dict <- FEATURE_OBJECT$get_scores_dict(matrix_dict$plddt, list_sequence_info)
    
    contact_matrix <- matrix_dict$contact_matrix
    
    # FEATURE_OBJECT.print_matrix_dict(matrix_dict)
    
    dc_obj <- domain_clustering(matrix_dict,
                                list_sequence_info,
                                alphafold_version = mode,
                                outdir = outdir, 
                                plotting = TRUE)
    dc_result <- dc_obj$run_domain_clustering()
    
    coevolutionary_domains <- dc_result$coevolutionary_domains
    interacting_coevolutionary_domains <- dc_result$interacting_coevolutionary_domains
    entity_region_dict <- dc_result$entity_region_dict
    
    INTERFACE_IDENTIFICATION <- interface_identification(interacting_coevolutionary_domains, 
                                                         entity_region_dict,
                                                         plddt_dict,
                                                         rec_sequence_list,
                                                         list_sequence_info,
                                                         contact_matrix,
                                                         contact_threshold,
                                                         chain_dict,
                                                         polymer_chain_dict)
    
    temp_interface <- INTERFACE_IDENTIFICATION$extract_interface()
    interface_dict <- temp_interface$interface_dict
    protein_interface_dict <- temp_interface$protein_interface_dict
    interaction_link_dict <- temp_interface$interaction_link_dict   
    
    temp_dfs <- INTERFACE_IDENTIFICATION$get_interface_info_dataframes(interface_dict, interaction_link_dict)
    interface_df_per_token <- temp_dfs$interface_df_per_token
    interface_df <- temp_dfs$interface_df
    
    ribbon_diagram <- RIBBON_DIAGRAM(list_sequence_info,
                                     interface_dict,
                                     protein_interface_dict,
                                     plddt_dict,
                                     outdir = outdir)
    
    ribbon_diagram$plot_ribbon_diagram()
    
    write_dataframe(interface_df_per_token, "interface_df_per_token", outdir)
    write_dataframe(interface_df, "binding_interfaces", outdir)
    
    return(list(interface_df_per_token = interface_df_per_token, interface_df = interface_df))
}

main <- function() {
    args <- parse_args()
    
    in_dir <- args$in_dir
    mode <- args$mode
    contact_threshold <- args$contact_threshold
    
    temp <- define_interfaces(in_dir, mode, contact_threshold)
    interface_df_per_token <- temp$interface_df_per_token
    interface_df <- temp$interface_df
    
    print("finished")
}

# Equivalent of Python's "if __name__ == '__main__':"
if (!interactive()) {
    main()
}


# Dependency Adaptation:
    # The Python "import" statements have been translated to either "library()" calls (for standard R packages) or "source()" calls for the external modules. Each "source()" is assumed to load an R script providing the equivalent functionality of the original Python module.
# Working Directory:
    # Python's "os.path.dirname(os.path.realpath(file))" is replaced with "getwd()" in R since file is not directly available. This approximates the current working directory.
# Argument Parsing:
    # The Python "argparse" module is replaced with the R package "optparse". We replicate the command-line argument behavior including printing help when no arguments are provided. The choices constraint for the "mode" argument is not enforced by optparse directly; it would need additional checking if necessary.
    # The "restricted_float" function is translated to an R function that performs the same type conversion and range checking.
# DataFrame Writing:
    # The "write_dataframe" function uses R's "write.csv" to write the dataframe to a CSV file with the same naming convention used in Python.
# Object Method Calls and Data Structures:
    # Python method calls such as "FEATURE_OBJECT.extract_feature_filepath()" are transposed to R using the "$" operator (e.g., FEATURE_OBJECT$extract_feature_filepath()). It is assumed that the sourced modules implement these as R objects with methods that return lists mimicking Python tuple returns.
# Domain Clustering and Interface Identification:
    # The results from functions like "domain_clustering" and "interface_identification" are expected to be returned as lists, with named elements corresponding to those in the Python version. Each function call that returns multiple values is unpacked accordingly.
# Script Execution:
    # The Python block checking "name == 'main'" is replaced by an R check that runs main() when the script is not running in interactive mode.
# Formatting:
    # All comments, indentation, and formatting have been preserved as closely as possible to the original Python code, with appropriate R syntax adjustments.
# Potential Issues:
    # The functionality of the external modules (e.g., CCM_AF3, compare_protein_seq, etc.) is assumed to have been appropriately translated to R and provided by the corresponding sourced scripts.
    # The optparse package in R does not enforce argument choices automatically (like the Python argparse "choices" parameter). Additional validation would be required if strict enforcement is needed.
