# ==============================
# Load Required Libraries
# ==============================
library(R6)
library(circlize)
library(RColorBrewer)
library(grDevices)
library(grid)

# Get the module directory (approximation of __file__ in R)
module_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% "module_placeholder.R"))

# ==============================
# Helper Functions
# ==============================
# Assign colours to interfaces
get_interface2colour <- function(interfaces_list) {
  palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(interfaces_list))
  return(setNames(colour_palette, interfaces_list))
}

# Assign colours to pLDDT values
get_colour_plddt <- function(plddt_value) {
  colour_breaks <- c(0, 50, 70, 90, 100)
  colours <- c("#ff7d45", "#ffdb13", "#65cbf3", "#0053d6")
  return(colours[cut(plddt_value, breaks = colour_breaks, labels = FALSE, include.lowest = TRUE)])
}

# Assign colours to conservation values
get_conservation_colour <- function(value) {
  colour_palette <- rev(brewer.pal(11, "RdBu"))
  return(colour_palette[round(value * (length(colour_palette) - 1)) + 1])
}

# Utility infix operator for default value
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# ==============================
# Ribbon Diagram Class
# ==============================
RIBBON_DIAGRAM <- R6Class("RIBBON_DIAGRAM",
  public = list(
    list_sequence_info = NULL,
    interface_dict = NULL,
    protein_interface_dict = NULL,
    plddt_dict = NULL,
    conservation_dict = NULL,
    outdir = NULL,

    # Initialise ribbon diagram class with sequence and interface data
    initialise <- function(list_sequence_info, interface_dict, protein_interface_dict,
                          plddt_dict = list(), conservation_dict = list(), outdir = "") {
      self$list_sequence_info <- list_sequence_info
      self$interface_dict <- interface_dict
      self$protein_interface_dict <- protein_interface_dict
      self$plddt_dict <- plddt_dict
      self$conservation_dict <- conservation_dict
      self$outdir <- outdir
    },

# Create ribbon diagram
create_ribbon_plot = function() {
  # Extract sequence information
  list_fasta_name    <- self$list_sequence_info[[1]]
  list_fasta_acclen  <- self$list_sequence_info[[2]]
  list_fasta_len     <- self$list_sequence_info[[4]]

  # Precompute interface colours
  cmap_plddt <- brewer.pal(11, "Spectral")
  cmap_conservation <- rev(brewer.pal(11, "RdBu"))

  get_plddt_colour_func <- function(value) cmap_plddt[round((value / 100) * (length(cmap_plddt) - 1)) + 1]
  get_conservation_colour <- function(value) cmap_conservation[round(value * (length(cmap_conservation) - 1)) + 1]

  interface2colour <- get_interface2colour(names(self$interface_dict))

  # Initialise Circos plot
  sectors <- setNames(as.numeric(list_fasta_len), list_fasta_name)
  circos.clear()
  circos.par(gap.degree = 5)
  circos.initialise(factors = names(sectors), xlim = matrix(c(rep(0, length(sectors)), sectors), ncol = 2))

  # Define sector sizes
  sectors <- setNames(as.numeric(list_fasta_len), list_fasta_name)
  circos.initialize(factors = names(sectors), xlim = matrix(c(rep(0, length(sectors)), sectors), ncol = 2))

  # Process each sector
  lapply(names(sectors), function(sector_name) {
  sector_size <- sectors[[sector_name]]

    # Fetch sequence-based values
    plddt_list <- self$plddt_dict[[sector_name]]
    conservation_list <- self$conservation_dict[[sector_name]]

    # Define track layout
    tracks_position_list <- if (length(conservation_list) > 0) list(c(75,85), c(88,93), c(95,100)) else list(c(75,85), c(88,93))

    # Add tracks
    lapply(tracks_position_list, function(tp) {
      circos.trackPlotRegion(factors = sector_name, ylim = c(0,1), track.height = 0.05,
                             panel.fun = function(x, y) circos.axis(h = "top", labels.cex = 0.6))
    })

    # Add pLDDT and conservation colour bars (vectorized)
    if (!is.null(plddt_list)) {
      circos.rect(xleft = 0:(sector_size - 1),
        ybottom = 0,
        xright = 1:sector_size,
        ytop = 1,
        sector.index = sector_name,
        track.index = 2,
        col = sapply(plddt_list, get_colour_plddt),
        border = NA)
    }

    if (!is.null(conservation_list)) {
      circos.rect(xleft = 0:(sector_size - 1),
        ybottom = 0,
        xright = 1:sector_size,
        ytop = 1,
        sector.index = sector_name,
        track.index = 3,
        col = sapply(conservation_list, get_conservation_colour),
        border = NA)
    }

    # Process interface ranges
    if (sector_name %in% names(self$protein_interface_dict)) {
      lapply(names(self$protein_interface_dict[[sector_name]]), function(interface_name) {
        interface_ranges <- self$protein_interface_dict[[sector_name]][[interface_name]]$interface_range
        col <- interface2colour[[interface_name]]
        lapply(interface_ranges, function(range) {
          circos.rect(xleft = range[1] - 1, ybottom = 0, xright = range[length(range)] - 1, ytop = 1,
                      sector.index = sector_name, track.index = 1, col = col, border = "black", lwd = 0.5)
        })
      })
    }
  })

  # Draw interface links (Vectorized)
  invisible(lapply(names(self$interface_dict), function(interface) {
    prot_1 <- self$interface_dict[[interface]]$prot_1$accesion_id
    prot_2 <- self$interface_dict[[interface]]$prot_2$accesion_id
    col <- adjustcolour(interface2colour[[interface]], alpha.f = 0.25)
    links <- self$interface_dict[[interface]]$links

    lapply(links, function(link) {
      circos.link(prot_1, c(link[[1]][1] - 1, link[[1]][2] - 1), prot_2, c(link[[2]][1] - 1, link[[2]][2] - 1), col = col, border = NA)
    })
  }))
}

# Save plot as PNG file
plot_ribbon_diagram = function() {
  circos_obj <- self$create_ribbon_plot()
  outdir <- self$outdir
  filename <- paste0(outdir, "/ribbon_plot.png")

  # Open a PNG device to save the plot
  png(filename, width = 1000, height = 1000, res = 300)
  # Re-draw the circos plot
  self$create_ribbon_plot()
  # Maybe a TIFF could be produced instead for lossless images?

  # Create legend using base graphics; position approximated to the right of the plot
  legend("topright", legend = c("Very high", "High", "Low", "Very Low"),
          fill = c("#0053d6", "#65cbf3", "#ffdb13", "#ff7d45"),
          title = "Model confidence",
          cex = 0.8,
          bty = "n")
  dev.off()
}
)
)
