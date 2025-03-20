# ==============================
# Load Required Libraries
# ==============================
library(optparse)   # For command-line argument parsing
library(tools)  # For file manipulation
library(data.table) # For efficient CSV writing with fwrite()

# Load external modules containing helper functions
# CCM_AF3; compare_protein_seq; domain_clustering; MMCIFPARSER, HSSPPARSER, alphafold_msa; conservation_score; interface_identification; RIBBON_DIAGRAM
modules <- c("confidence_contact_matrix.R",
            "alignment_utils.R",
            "domain_clustering.R", 
             "parsers.R",
             "interface_identification.R",
             "ribbon_diagram.R")
for (module in modules) {
    source(file.path("src/module", module))
}

# Set working directory
working_dir <- getwd()

# ==============================
# Command-Line Argument Parsing
# ==============================
parse_args <- function() {
    # Define available command-line arguments
    option_list <- list(
        make_option(c("-i"), dest = "in_dir", type = "character", default = "",
                    help = "Path to a directory where input folder is stored."),
        make_option(c("-c"), dest = "config_dir", type = "character", default = "",
                    help = "Path to a directory where configuration files are stored."),
        make_option(c("-m", "--mode"), dest = "mode", type = "character", default = "AF3",
                    help = "AlphaFold output type. Options: AF3, AF2, ColabFold."),
        make_option(c("-t", "--threshold"), dest = "contact_threshold", type = "double", default = 0.7,
                    help = "Threshold for detecting a contact-link in the contact probability matrix (range: 0.0 - 1.0).")
    )

    parser <- OptionParser(option_list = option_list)

    # Ensure user provides arguments
    args_cmd <- commandArgs(trailingOnly = TRUE)
    if (length(args_cmd) == 0) {
        print_help(parser)
        stop("No arguments provided.")
    }

    args <- parse_args(parser)

    # Validate mode argument
    if (!(args$mode %in% c("AF3", "AF2", "ColabFold"))) {
        stop("Invalid mode. Choose from: AF3, AF2, ColabFold.")
    }

    # Ensure threshold is valid
    if (is.null(args$contact_threshold)) args$contact_threshold <- 0.7
    args$contact_threshold <- restricted_float(args$contact_threshold)

    return(args)
}

# ==============================
# Helper Functions
# ==============================
# Validate floating-point values within range [0.0, 1.0]
restricted_float <- function(x) {
    # Try to convert x to numeric
    x_numeric <- as.numeric(x)
    if (is.na(x_numeric)) stop(sprintf("'%s' not a floating-point literal", x))
    if (x_numeric < 0.0 || x_numeric > 1.0) stop(sprintf("'%s' not in range [0.0, 1.0]", x))
    return(x_numeric)
}

#Save DataFrame to CSV with error handling
write_dataframe <- function(df, filename, outdir_path) {
    filepath <- file.path(outdir_path, paste0(filename, ".csv"))
        tryCatch({
            fwrite(df, filepath)
        }, error = function(e) {
            message("Error writing file:", filepath)
            stop(e)
        })
}

# ==============================
# Interface Detection
# ==============================
define_interfaces <- function(in_dir, mode, contact_threshold) {
    outdir <- file.path(in_dir, "AlphaBridge")

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    # Feature extraction for AF3
    if (mode == "AF3") {
        FEATURE_OBJECT <- CCM_AF3(in_dir)

        list(feature_path, structure_path, job_request_path) <- FEATURE_OBJECT$extract_feature_filepath()

        list(list_sequence_info, rec_sequence_list, structure_sequence_list, polymer_chain_dict) <- FEATURE_OBJECT$extract_sequence_info()

    } else {
        stop("Output from AF2 or ColabFold not implemented yet")
    }

    # Compare sequences and extract chain information
    chain_dict <- compare_protein_seq(structure_sequence_list, rec_sequence_list)$extract_chain_dict()

    # Extract confidence scores
    matrix_dict <- FEATURE_OBJECT$extract_matrix_dict()
    plddt_dict <- FEATURE_OBJECT$get_scores_dict(matrix_dict$plddt, list_sequence_info)
    contact_matrix <- matrix_dict$contact_matrix

    # FEATURE_OBJECT.print_matrix_dict(matrix_dict)

    dc_obj <- domain_clustering(matrix_dict,
                                list_sequence_info,
                                alphafold_version = mode,
                                outdir = outdir, 
                                plotting = TRUE)
    # Perform domain clustering
    dc_result <- dc_obj$run_domain_clustering()

    coevolutionary_domains <- dc_result$coevolutionary_domains
    interacting_coevolutionary_domains <- dc_result$interacting_coevolutionary_domains
    entity_region_dict <- dc_result$entity_region_dict

    # Identify interfaces
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

    # Generate ribbon diagrams
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

# Run main function
if (!interactive()) {
    main()
}
