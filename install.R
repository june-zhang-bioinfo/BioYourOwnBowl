# install.R
# Installation script for BioYourOwnBowl

message("========================================")
message("Installing BioYourOwnBowl")
message("========================================\n")

# Setup library
sc_lib <- "~/R/bioyourownbowl_library"

if (dir.exists(sc_lib)) {
  message("Library already exists at: ", sc_lib)
  response <- readline(prompt = "Reinstall from scratch? This will delete existing packages. (y/N): ")
  if (tolower(response) == "y") {
    message("Removing existing library...")
    unlink(sc_lib, recursive = TRUE)
    dir.create(sc_lib, recursive = TRUE)
  } else {
    message("Using existing library. Will update/install missing packages only.")
  }
} else {
  message("Creating new library at: ", sc_lib)
  dir.create(sc_lib, recursive = TRUE)
}

# Set library paths
.libPaths(c(sc_lib, .libPaths()))

message("\nLibrary paths:")
message("  [1] ", .libPaths()[1], " (BioYourOwnBowl)")
message("  [2] ", .libPaths()[2], " (shared dependencies)")

# Check for remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
  message("\nInstalling remotes...")
  install.packages("remotes")
}

# Install dependencies
message("\n--- Installing GitHub dependencies ---")

install_gh <- function(repo, lib) {
  pkg_name <- sub(".*/", "", sub("@.*", "", repo))
  message("Installing ", pkg_name, "...")
  tryCatch({
    remotes::install_github(repo, lib = lib, upgrade = "never")
    message("  ✓ ", pkg_name, " installed")
  }, error = function(e) {
    message("  ✗ ", pkg_name, " failed: ", e$message)
  })
}

install_gh("satijalab/seurat@v5.1.0", sc_lib)
install_gh("tidyverse/ggplot2@v3.5.2", sc_lib)
install_gh("settylab/convert2anndata", sc_lib)
install_gh("jokergoo/ComplexHeatmap", sc_lib)
install_gh("kevinblighe/EnhancedVolcano", sc_lib)

message("\n--- Installing BioYourOwnBowl ---")
install_gh("june-zhang-bioinfo/BioYourOwnBowl", sc_lib)

# Verify installation
message("\n--- Verifying installation ---")
if (requireNamespace("BioYourOwnBowl", quietly = TRUE)) {
  message("✓ BioYourOwnBowl version ", packageVersion("BioYourOwnBowl"))
  message("✓ Seurat version ", packageVersion("Seurat"))
  message("✓ ggplot2 version ", packageVersion("ggplot2"))
  
  message("\n========================================")
  message("Installation complete!")
  message("========================================\n")
  message("To use BioYourOwnBowl, start your script with:")
  message('.libPaths(c("~/R/bioyourownbowl_library", .libPaths()))')
  message('library(BioYourOwnBowl)')
} else {
  message("\n✗ Installation may have failed. Check errors above.")
}