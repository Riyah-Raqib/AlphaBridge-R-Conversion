# ==============================
# Load Required Libraries
# ==============================
library(httr)
library(jsonlite)
library(R6)

# Attempt to load 'stringr' for regex operations (using base R regex if unavailable)
library(stringr)

# Determine the module directory.
module_dir <- tryCatch({
  # In R, there is no direct equivalent of __file__, so we attempt to mimic its behavior.
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) {
  "."
})

WEBSITE_API <- "https://alphafold.ebi.ac.uk/api/prediction"
KEY_API <- "key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94"

ALPHAMISSENSE <- R6Class("ALPHAMISSENSE",
  public = list(
    UNIPROT_ID = NULL,

    initialise = function(UNIPROT_ID) {
      self$UNIPROT_ID <- UNIPROT_ID
    },

    extract_am_annotation = function() {
      UNIPROT_ID <- self$UNIPROT_ID

      r <- get_url(paste0(WEBSITE_API, "/", UNIPROT_ID, "?", KEY_API))

      # Parse JSON and take the first element
      data <- fromJSON(content(r, "text", encoding = "UTF-8"))[[1]]

      am_annotation <- read.csv(data$amAnnotationsUrl, stringsAsFactors = FALSE)

      return(am_annotation)
    },

    get_pathogenicity_list = function() {
      am_annotation <- self$extract_am_annotation()

      # Apply the split_column function to each row of the data frame
      am_annotation_modified <- do.call(rbind, lapply(1:nrow(am_annotation), function(i) {
        # Convert row to list for processing
        row <- as.list(am_annotation[i, , drop = FALSE])
        split_column(row)
      }))

      # Convert factors to characters if necessary
      if (is.factor(am_annotation_modified$REF)) {
        am_annotation_modified$REF <- as.character(am_annotation_modified$REF)
      }
      if (is.factor(am_annotation_modified$POS)) {
        am_annotation_modified$POS <- as.character(am_annotation_modified$POS)
      }

      # Group by 'REF' and 'POS' and calculate mean of 'am_pathogenicity'
      # Since there is no direct groupby command in base R, we use aggregate.
      grouped <- aggregate(am_pathogenicity ~ REF + POS, data = am_annotation_modified, FUN = mean)

      # Convert POS to numeric for proper sorting, mimicking int() conversion in Python
      grouped$numeric_POS <- as.numeric(grouped$POS)

      # Sort by numeric_POS
      grouped_sorted <- grouped[order(grouped$numeric_POS), ]

      # Extract the means as the pathogenicity list
      pathogenicity_list <- grouped_sorted$am_pathogenicity

      return(pathogenicity_list)
    }
  )
)

get_url <- function(url, ...) {
  response <- GET(url, ...)

  if (!(response$status_code >= 200 && response$status_code < 300)) {
    # Print the response text if the request was not successful
    print(content(response, "text", encoding = "UTF-8"))
    stop("HTTP error: ", response$status_code)
  }

  return(response)
}

split_column <- function(row) {
  protein_variant <- row$protein_variant

  # Use regex to find all matches of letters or digits
  # The regex pattern "[A-Za-z]+|\\d+" finds one or more letters or a digit sequence
  matches <- regmatches(protein_variant, gregexpr("[A-Za-z]+|\\d+", protein_variant))[[1]]

  # Assign the extracted values to new keys 'REF', 'POS', and 'ALT'
  row$REF <- matches[1]
  row$POS <- matches[2]
  row$ALT <- matches[3]

  return(row)
}
