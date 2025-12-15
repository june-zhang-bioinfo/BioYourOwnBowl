#' Plot Features with Dynamic Layout and Cutoffs
#'
#' Generates FeaturePlots for selected features from a Seurat object, 
#' adjusting the plot grid dynamically based on the number of features.
#' Saves both a scaled and an original version of the plot as PNGs.
#'
#' @param object A Seurat object.
#' @param features Character vector of features to plot.
#' @param n_col Integer. Number of columns in the plot grid. Default is 4.
#' @param pt.size Point size. Default is 0.8.
#' @param base_width Width of each plot in pixels. Default is 1000.
#' @param base_height Height of each plot in pixels. Default is 1000.
#' @param percentile The percentile (between 0 and 1) of feature values to be shown. The higher the percentile, the sharper the feature. Default is 0.9.
#' @param file_name Character string for output png name. If NULL, plot will be printed in console but not saved.
#'   Default is NULL.
#'
#' @return NULL. Plots are saved as PNG files in the working directory.
#' @export
features_plots <- function(object,
                           features,
                           n_col = 4,
                           pt.size = 0.8,
                           base_width = 1000,
                           base_height = 1000,
                           percentile = 0.9,
                           file_name = NULL
                           ) {
  
  # Ensure selected features exist in the object
  features_to_plot <- intersect(features, c(rownames(object), colnames(object@meta.data)))
  n_features <- length(features_to_plot)
  if (n_features == 0) stop("No selected features found in the Seurat object.")
  
  p_original <- FeaturePlot_scCustom(object, features = features_to_plot, num_columns = n_col, pt.size = pt.size)
  
  # Calculate grid layout
  n_row <- ceiling(n_features / n_col)
  width_scaled <- base_width * n_col
  height_scaled <- base_height * n_row
  
  # Save original plot
  if (!is.null(file_name)) {
    png(paste0(file_name, "_original.png"), width = width_scaled, height = height_scaled, res = 300)
    print(p_original)
    dev.off()
  }
  
  # Calculate cutoffs per feature (user-specified percentile)
  percentiles <- sapply(features_to_plot, function(g) {
    q <- quantile(FetchData(object, vars = g), percentile, na.rm = TRUE, layer = "data")
    if (q == 0) q <- 1e-9
    return(q)
  })
  
  # Create individual scaled plots with custom cutoff
  plots <- purrr::map2(
    features_to_plot,
    percentiles,
    ~ FeaturePlot_scCustom(
      object,
      pt.size = pt.size,
      features = .x,
      na_cutoff = .y,
      order = TRUE,
      raster = FALSE
    )
  )
  
  p_scaled_wrapped <- patchwork::wrap_plots(plots, ncol = n_col)
  
  # Save scaled plots
  if (!is.null(file_name)) {
    png(paste0(file_name, "_scaled.png"), width = width_scaled, height = height_scaled, res = 300)
    print(p_scaled_wrapped)
    dev.off()
  }
  
  return(list(
    original = p_original,
    scaled = p_scaled_wrapped
  ))
}


#' Create Violin Plot with Boxplot and Dotplot Overlay
#'
#' Generate a combined violin, boxplot, and dotplot visualization for gene expression
#' across different cell identities in a Seurat object with optional statistical comparisons.
#'
#' @param object A Seurat object containing expression data
#' @param features Character vector of gene names to plot
#' @param idents Character string specifying the metadata column to use for grouping cells.
#'   This will be set as the active identity of the Seurat object.
#' @param groups Character vector specifying which identity groups to include and their order.
#'   If NULL, all identities will be used. Default is NULL.
#' @param colors Named vector of colors for each identity group. If NULL, default ggplot2 colors
#'   will be used. Default is NULL.
#' @param binwidth Numeric value for dotplot bin width. Default is 0.2.
#' @param file_name Character string for output png name. If NULL, plot will be displayed but not saved.
#'   Default is NULL.
#' @param width Numeric value for plot width in inches. Default is 8.
#' @param height Numeric value for plot height in inches. If NULL, calculated dynamically.
#' @param ncol Integer specifying number of columns in facet layout. Default is 5.
#' @param add_stats Logical indicating whether to add statistical test results to the plot.
#'   Default is FALSE.
#' @param test_method Character string specifying the statistical test method. Options are
#'   "wilcox" (Wilcoxon rank-sum test) or "t.test" (t-test). Default is "wilcox".
#' @param comparisons List of length-2 character vectors specifying pairwise comparisons to test.
#'   If NULL and add_stats = TRUE, all pairwise comparisons will be performed. Default is NULL.
#' @param p_adjust_method Character string specifying p-value adjustment method for multiple
#'   testing. Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'   Default is "BH" (Benjamini-Hochberg).
#' @param symnum_args List of arguments to pass to symnum() for significance symbols.
#'   Default is list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
#'                   symbols = c("****", "***", "**", "*", "ns")).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Basic usage without statistics
#' violin_plots(
#'   object = object,
#'   features = c("CD3D", "CD8A", "CD4"),
#'   idents = "cell_type"
#' )
#'
#' # With statistics - all pairwise comparisons
#' violin_plots(
#'   object = object,
#'   features = gene_list,
#'   idents = "cohort_tp",
#'   groups = c("Chronic Post-DAA", "ACTG Post-DAA", "SR ≥24W-post-resolution"),
#'   add_stats = TRUE
#' )
#'
#' # With specific comparisons
#' my_comparisons <- list(
#'   c("Chronic Post-DAA", "ACTG Post-DAA"),
#'   c("Chronic Post-DAA", "SR ≥24W-post-resolution")
#' )
#'
#' violin_plots(
#'   object = object,
#'   features = gene_list,
#'   idents = "cohort_tp",
#'   groups = c("Chronic Post-DAA", "ACTG Post-DAA", "SR ≥24W-post-resolution"),
#'   colors = c("Chronic Post-DAA" = "#F8A19FB2",
#'              "ACTG Post-DAA" = "#2ED9FFB2",
#'              "SR ≥24W-post-resolution" = "#FEAF16B2"),
#'   add_stats = TRUE,
#'   comparisons = my_comparisons,
#'   test_method = "wilcox",
#'   file_name = "expression_violin.png"
#' )
#' }
#'
#' @export
violin_plots <- function(object,
                         features,
                         idents,
                         groups = NULL,
                         colors = NULL,
                         binwidth = 0.2,
                         file_name = NULL,
                         width = 8,
                         height = NULL,
                         ncol = 5,
                         add_stats = FALSE,
                         test_method = "wilcox",
                         comparisons = NULL,
                         p_adjust_method = "BH",
                         symnum_args = list(
                           cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                           symbols = c("****", "***", "**", "*", "ns")
                         )) {
  
  # Set active identity
  levels <- object[[idents]] %>% unlist() %>% as.character() %>% stringr::str_sort() %>% unique()
  Idents(object) <- object[[idents]] %>% unlist() %>% as.character()
  Idents(object) <- factor(object[[idents]] %>% unlist() %>% as.character(), levels = levels)
  
  # Subset cells if groups specified
  if (!is.null(groups)) {
    cells_to_use <- Seurat::WhichCells(object, idents = groups)
  } else {
    cells_to_use <- colnames(object)
    groups <- levels(Idents(object))
  }
  
  # Fetch expression data
  features <- intersect(features, c(rownames(object), colnames(object@meta.data)))
  expr_mat <- Seurat::FetchData(object, vars = features, cells = cells_to_use, layer = "data")
  expr_mat$ident <- Seurat::Idents(object)[cells_to_use]
  
  # Convert to long format
  expr_long <- expr_mat %>%
    rownames_to_column(var = "cell") %>%
    pivot_longer(cols = all_of(features), names_to = "gene", values_to = "expression")
  
  # Set factor levels
  expr_long$ident <- factor(expr_long$ident, levels = groups)
  expr_long$gene <- factor(expr_long$gene, levels = unique(expr_long$gene))
  
  # Calculate dynamic height if not provided
  if (is.null(height)) {
    n_genes <- length(features)
    n_rows <- ceiling(n_genes / ncol)
    # Add extra height if statistics are being added
    height_per_row <- ifelse(add_stats, 3.5, 3)
    height <- max(3, n_rows * height_per_row)
  }
  
  # Generate all pairwise comparisons if not provided
  if (add_stats && is.null(comparisons)) {
    if (length(groups) > 1) {
      comparisons <- utils::combn(groups, 2, simplify = FALSE)
    } else {
      warning("Cannot perform statistical tests with only one group")
      add_stats <- FALSE
    }
  }
  
  # Create plot
  p <- ggplot(expr_long, aes(x = ident, y = expression, color = ident)) +
    geom_violin(
      data = subset(expr_long, expression >= 0),
      trim = TRUE,
      scale = "count",
      adjust = 1,
      alpha = 0.7
    ) +
    geom_dotplot(
      data = expr_long,
      binaxis = 'y',
      stackdir = 'center',
      dotsize = 0.35,
      binwidth = binwidth,
      method = "histodot"
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.3,
      color = "black",
      fill = "white"
    ) +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      strip.text = element_text(size = 9),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Identity", y = "Expression")
  
  # Add statistical comparisons
  if (add_stats) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = test_method,
      label = "p.signif",
      symnum.args = symnum_args,
      p.adjust.method = p_adjust_method
    )
  }
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  # Save plot if filename provided
  if (!is.null(file_name)) {
    ggsave(file_name, plot = p, height = height, width = width)
  }
  
  return(p)
}



#' Save Split DimPlots from a Seurat Object
#'
#' This function generates and saves UMAP/tSNE DimPlots from a Seurat object, 
#' split by a specified metadata variable. The plots are arranged automatically 
#' based on the number of levels in the splitting variable.
#'
#' @param object A \code{Seurat} object.
#' @param split_var A character string specifying the metadata column to split plots by.
#' @param group_var A character string specifying the metadata column to color (group) cells by.
#' @param colors A character vector of colors to use for the groups.
#' @param file_name Character string for output png name. If NULL, plot will be displayed but not saved.
#'   Default is NULL.
#'
#' @return Invisibly returns the file path of the saved plot.
#' @export
#'
#' @examples
#' \dontrun{
#' batch_vars <- "Cohort"
#' my_colors <- c("red", "blue", "green")
#' dimplots_split(
#'   object = object,
#'   split_var = batch_vars,
#'   group_var = "seurat_clusters",
#'   colors = my_colors,
#'   file_name = NULL
#' )
#' }
#'
dimplots_split <- function(object, split_var, group_var, colors, file_name = NULL) {
  if (length(split_var) != 1) {
    stop("split_var must be a single metadata column name (character string).")
  }
  
  num_levels <- length(unique(object@meta.data[[split_var]]))
  
  num_columns <- ceiling(sqrt(num_levels))
  num_rows <- ceiling(num_levels / num_columns)
  
  width <- num_columns * 3
  height <- num_rows * 3
  
  p <- Seurat::DimPlot(
    object,
    split.by = split_var,
    group.by = group_var,
    label = FALSE,
    ncol = num_columns,
    cols = colors
  )
  
  if(!is.null(file_name)){
    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
  }

  return(p)
}


#' Set up Python environment for BioYourOwnBowl
#' @export
setup_python <- function() {
  # Check for Python with numpy
  cfg <- reticulate::py_discover_config(required_module = "numpy", use_environment = NULL)
  
  # Install Python if none found OR if numpy isn't available
  if (is.null(cfg$python) || is.null(cfg$numpy)) {
    message("No suitable Python found. Installing Python 3.10 via reticulate...")
    reticulate::install_python(version = "3.10.11")
    # Update config after installation
    cfg <- reticulate::py_discover_config()
  }
  
  env_name <- "r-reticulate"
  env_path <- file.path(reticulate::virtualenv_root(), env_name)
  
  if (!dir.exists(env_path)) {
    message("Creating virtual environment r-reticulate...")
    reticulate::virtualenv_create(env_name, python = cfg$python)
    message("Installing packages: scanpy, anndata...")
    reticulate::virtualenv_install(env_name, c("scanpy", "anndata"))
  }
  
  reticulate::use_virtualenv(env_name, required = TRUE)
  
  # Verify installation
  tryCatch({
    reticulate::import("scanpy")
    message("Python environment setup complete!")
  }, error = function(e) {
    warning("Setup completed but scanpy import failed: ", e$message)
  })
}



#' Prepare a Seurat Object for Export to AnnData (.h5ad)
#'
#' This function subsets, normalizes, and converts a Seurat object into
#' an AnnData-compatible format for downstream analysis in Python (e.g. Scanpy).
#' It optionally filters metadata and writes the result to an `.h5ad` file.
#'
#' @param object A Seurat object.
#' @param features Character vector of feature (gene) names to retain.
#' @param metadata_vars Character vector of metadata column names to keep.
#' @param file_name Character string giving the output `.h5ad` file name. Default is anndata.
#' @param assay Character string giving the assay to use. Defaults to
#'   the Seurat object's current default assay.
#' @param normalize Logical; whether to run \code{Seurat::NormalizeData()}.
#'   Defaults to \code{TRUE}.
#'
#' @return Invisibly returns the output file path.
#' @export
#'
#' @examples
#' \dontrun{
#' features <- unique(c(VariableFeatures(object), features))
#' metadata_vars <- c("cell_id", names(density_plot_specs))
#' prepare_h5ad(object, features, metadata_vars, "exported_data.h5ad")
#' }
prepare_h5ad <- function(
    object,
    features,
    metadata_vars,
    file_name = "anndata",
    assay = NULL,
    normalize = TRUE
) {
  cnmf <- subset(object, features = features)
  
  if (normalize) {
    cnmf <- Seurat::NormalizeData(cnmf)
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(cnmf)
  }
  
  cnmf@assays[[assay]]$counts <- cnmf@assays[[assay]]$data
  cnmf@assays[[assay]]$data <- NULL
  cnmf@assays[[assay]]$scale.data <- NULL
  
  cnmf@meta.data <- subset(cnmf@meta.data, select = metadata_vars)
  
  sce <- convert_seurat_to_sce(cnmf)
  ad <- convert_to_anndata(sce, assayName = "data", useAltExp = TRUE)
  
  write_h5ad(ad, file_name)
  
  message("Saved AnnData file: ", file_name)
  invisible(file_name)
}



#' Generate Python Code for Scanpy UMAP Density Plots
#'
#' This function generates a block of Python code that, when executed,
#' produces UMAP density plots for specified metadata variables in an
#' AnnData object. It is designed for cross-language Seurat → Scanpy
#' workflows, where an R object is exported to `.h5ad` format and
#' processed in Python using `scanpy`.
#'
#' @param file_name Character string. Path to the `.h5ad` file that contains
#'   the AnnData object to be plotted.
#' @param density_specs A named list specifying plotting order for each metadata
#'   variable. Each element should itself be a list containing an `order`
#'   vector (e.g., `list(order = c("group1", "group2", ...))`).
#' @param base_figsize Numeric vector of length 2 giving the base figure size
#'   in inches for each density plot panel. Default is `c(4, 4)`.
#' @param max_per_row Integer specifying the maximum number of plots per row.
#'   Used to determine the number of columns in multi-panel figures.
#'
#' @return A single character string containing executable Python code.
#'   This code can be passed directly to `reticulate::py_run_string()`.
#'
#' @details
#' The generated Python code will:
#' \enumerate{
#'   \item Import required libraries (`scanpy` and `matplotlib`).
#'   \item Read the `.h5ad` file into an AnnData object.
#'   \item For each metadata variable in \code{density_specs}, compute and plot
#'         UMAP density maps using `scanpy.tl.embedding_density()` and
#'         `scanpy.pl.embedding_density()`.
#'   \item Save each density plot to a PNG file named
#'         \code{density_<meta_col>.png}.
#' }
#'
#' @examples
#' \dontrun{
#' density_plot_specs <- list(
#'   cohort_tp = list(order = c("Chronic Pre-DAA", "Chronic Post-DAA")),
#'   protein_group = list(order = c("Group1", "Group2"))
#' )
#'
#' py_code <- density_plot(
#'   file_name = "mydata.h5ad",
#'   density_specs = density_plot_specs
#' )
#'
#' }
#'
#' @seealso \code{\link[reticulate]{py_run_string}}
#'
#' @export
density_plot <- function(file_name, density_specs,
                                       base_figsize = c(4, 4),
                                       max_per_row = 4) {
  code_blocks <- c(
    "import scanpy as sc",
    "import matplotlib.pyplot as plt",
    glue::glue("adata = sc.read_h5ad('{file_name}')")
  )
  
  for (meta_col in names(density_specs)) {
    groups <- density_specs[[meta_col]]$order
    # ncols <- ceiling(length(groups) / max_per_row)
    ncols <- ifelse(length(groups) < max_per_row, length(groups), max_per_row)
    
    fig_width <- base_figsize[1]
    fig_height <- base_figsize[2]
    
    code_block <- glue::glue("
sc.tl.embedding_density(adata, basis='umap', groupby='{meta_col}')
with plt.rc_context({{'figure.figsize': ({fig_width}, {fig_height})}}):
    sc.pl.embedding_density(
        adata,
        basis='umap',
        groupby='{meta_col}',
        bg_dotsize=20,
        color_map='inferno',
        fg_dotsize=80,
        ncols={ncols},
        group={jsonlite::toJSON(groups, auto_unbox=TRUE)},
        show=False
    )
plt.savefig('density_{meta_col}.png', bbox_inches='tight', dpi=300)
plt.close('all')
")
    code_blocks <- c(code_blocks, code_block)
  }
  
  reticulate::py_run_string(paste(code_blocks, collapse = "\n"))
  
  return(paste(code_blocks, collapse = "\n"))
}


#' Create Volcano Plot from Differential Expression Analysis
#'
#' This function performs differential expression analysis between two groups
#' and creates volcano plots with optional gene labeling using EnhancedVolcano.
#'
#' @param object A Seurat object containing expression data and metadata
#' @param idents Character string specifying the metadata column containing
#'   group identities. This column will be set as the active identity.
#' @param ident1 Character string specifying the first identity group (numerator
#'   in the comparison, typically the "test" group).
#' @param ident2 Character string specifying the second identity group (denominator
#'   in the comparison, typically the "control" group).
#' @param show_labels Can be one of three options:
#'   \itemize{
#'     \item FALSE (default): No gene labels are shown
#'     \item TRUE: All gene labels are shown
#'     \item Character vector: Only specified genes are labeled
#'   }
#' @param log2fc_cutoff Numeric value for log2 fold change threshold to highlight
#'   significant genes. Default is 1.
#' @param p_cutoff Numeric value for adjusted p-value cutoff to highlight significant
#'   genes. Default is 0.01.
#' @param file_name Optional character string for output file prefix. If provided,
#'   plots will be saved as PNG files. If NULL, plots are displayed but not saved.
#'   Default is NULL.
#' @param title Optional character string for plot title. If NULL, no title is shown.
#'   Default is NULL.
#' @param width Numeric value for plot width in pixels when saving. Default is 2000.
#' @param height Numeric value for plot height in pixels when saving. Default is 2000.
#' @param resolution Numeric value for plot resolution (DPI) when saving. Default is 300.
#' @param point_size Numeric value for point size in the volcano plot. Default is 1.
#' @param label_size Numeric value for gene label text size. Default is 2.
#' @param boxed_labels Logical indicating whether to draw boxes around gene labels.
#'   Default is TRUE.
#' @param draw_connectors Logical indicating whether to draw connectors from points
#'   to labels. Default is TRUE.
#' @param max_overlaps Integer specifying maximum number of overlapping labels allowed.
#'   Default is 50.
#'
#' @return A list containing:
#' \describe{
#'   \item{markers}{Data frame with differential expression results from FindMarkers}
#'   \item{n1}{Number of cells in ident1}
#'   \item{n2}{Number of cells in ident2}
#'   \item{deg_up}{Number of upregulated genes (log2FC > cutoff, p.adj < cutoff)}
#'   \item{deg_down}{Number of downregulated genes (log2FC < -cutoff, p.adj < cutoff)}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Counts cells in each group
#'   \item Performs differential expression using Seurat::FindMarkers
#'   \item Replaces p-values of 0 with machine epsilon to avoid plotting issues
#'   \item Creates volcano plot(s) with EnhancedVolcano
#'   \item Optionally saves plot(s) as PNG files
#' }
#'
#' @examples
#' \dontrun{
#' # No labels (default)
#' result <- plot_volcano_de(
#'   object = seurat_obj,
#'   idents = "cell_type",
#'   ident1 = "CD8 T",
#'   ident2 = "CD4 T"
#' )
#'
#' # Show all labels
#' result <- plot_volcano_de(
#'   object = seurat_obj,
#'   idents = "treatment",
#'   ident1 = "Treated",
#'   ident2 = "Control",
#'   show_labels = TRUE
#' )
#'
#' # Show specific genes only
#' selected <- c("GNLY", "GZMB", "PRF1", "CCL4")
#' result <- plot_volcano_de(
#'   object = seurat_obj,
#'   idents = "Timepoint disease stage",
#'   ident1 = "FC",
#'   ident2 = "Acute HBV resolved",
#'   show_labels = selected,
#'   log2fc_cutoff = 1,
#'   p_cutoff = 0.01,
#'   file_name = "volcano/FC_vs_Acute",
#'   title = "CD8 T cells"
#' )
#' }
#'
#' @import Seurat
#' @import EnhancedVolcano
#'
#' @export
volcano_plots <- function(object,
                          idents,
                          ident1,
                          ident2,
                          show_labels = FALSE,
                          log2fc_cutoff = 1,
                          p_cutoff = 0.01,
                          file_name = NULL,
                          title = NULL,
                          width = 2000,
                          height = 2000,
                          resolution = 300,
                          point_size = 1,
                          label_size = 2,
                          boxed_labels = TRUE,
                          draw_connectors = TRUE,
                          max_overlaps = 50) {
  
  # Validate inputs
  if (!idents %in% colnames(object@meta.data)) {
    stop("Column '", idents, "' not found in metadata")
  }
  
  # if(is.null(file_name)){file_name <- paste0(ident1, "_vs_", ident2)}
  
  # Count cells in each group
  n1 <- sum(object[[idents, drop = TRUE]] == ident1, na.rm = TRUE)
  n2 <- sum(object[[idents, drop = TRUE]] == ident2, na.rm = TRUE)
  
  if (n1 == 0 || n2 == 0) {
    stop("One or both groups have zero cells. Check ident1 and ident2 values.")
  }
  
  # Set active identity
  Idents(object) <- object[[idents]] %>% unlist() %>% as.character()
  
  # Perform differential expression
  message("Performing differential expression: ", ident1, " vs. ", ident2)
  markers <- FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = 0,
    min.pct = 0
  )
  
  markers$gene <- rownames(markers)
  markers$p_val[markers$p_val == 0] <- .Machine$double.xmin
  
  deg <- subset(markers, abs(avg_log2FC) > log2fc_cutoff & p_val_adj < p_cutoff)
  
  # Count DEGs
  deg_up <- sum(markers$avg_log2FC > log2fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE)
  deg_down <- sum(markers$avg_log2FC < -log2fc_cutoff & markers$p_val_adj < p_cutoff, na.rm = TRUE)
  
  # Prepare caption
  caption_text <- paste0(
    n1, " cells vs. ", n2, " cells\n",
    deg_up, " genes ↑ & ", deg_down, " genes ↓"
  )
  
  # Determine label settings based on show_labels
  if (is.logical(show_labels)) {
    if (show_labels) {
      # TRUE: show all labels
      lab_size <- label_size
      select_lab <- NULL
      boxed <- boxed_labels
      connectors <- draw_connectors
    } else {
      # FALSE: no labels
      lab_size <- 0
      select_lab <- NULL
      boxed <- FALSE
      connectors <- FALSE
    }
  } else if (is.character(show_labels)) {
    # Character vector: show specific genes
    lab_size <- label_size
    select_lab <- intersect(show_labels, deg$gene)
    boxed <- boxed_labels
    connectors <- draw_connectors
  } else {
    stop("show_labels must be TRUE, FALSE, or a character vector of gene names")
  }
  
  # Create volcano plot function
  create_volcano <- function(with_labels = TRUE) {
    if (!with_labels) {
      # No labels version
      EnhancedVolcano::EnhancedVolcano(
        markers,
        lab = rownames(markers),
        boxedLabels = FALSE,
        labSize = 0,
        pointSize = point_size,
        lengthConnectors = unit(0.01, "npc"),
        pCutoff = p_cutoff,
        x = "avg_log2FC",
        y = "p_val_adj",
        FCcutoff = log2fc_cutoff,
        title = title,
        subtitle = paste0(ident1, " vs. ", ident2),
        subtitleLabSize = 15,
        caption = caption_text,
        drawConnectors = FALSE,
        max.overlaps = max_overlaps
      )
    } else {
      # With labels version
      EnhancedVolcano::EnhancedVolcano(
        markers,
        lab = rownames(markers),
        boxedLabels = boxed,
        labSize = lab_size,
        pointSize = point_size,
        lengthConnectors = unit(0.01, "npc"),
        pCutoff = p_cutoff,
        x = "avg_log2FC",
        y = "p_val_adj",
        FCcutoff = log2fc_cutoff,
        title = title,
        subtitle = paste0(ident1, " vs. ", ident2),
        subtitleLabSize = 15,
        caption = caption_text,
        selectLab = select_lab,
        drawConnectors = connectors,
        max.overlaps = max_overlaps
      )
    }
  }
  
  if (is.character(show_labels)) {
    p <- create_volcano(with_labels = TRUE)
  } else {
    p <- create_volcano(with_labels = isTRUE(show_labels))
  }
  
  print(p)
  
  # If show_labels is a character vector, save both versions
  if(!is.null(file_name)){
    png(file_name, res = resolution, height = height, width = width)
    grid::grid.newpage()
    grid::grid.draw(p)
    dev.off()
    message("Saved: ", file_name)
  }

  return(markers)
}

#' Generate 2D Log2 Fold-Change Comparison Plots Between Identity Pairs
#'
#' This function performs pairwise differential expression analysis between
#' specified identity groups in a Seurat object, merges results from two
#' comparisons, classifies genes by directionality and significance, and
#' visualizes the results as 2D scatter plots of \eqn{log2FC_1} vs. \eqn{log2FC_2}.
#'
#' @param object A \code{Seurat} object containing expression and metadata.
#' @param comparison_table A data frame specifying pairs of identities to compare.
#'   Each row should include four columns (e.g., \code{ident1}, \code{ident2},
#'   \code{ident3}, \code{ident4}) defining two comparisons.
#' @param output_dir Character. Path to the directory where plots will be saved.
#'   Defaults to \code{"2dim_log2fc"}.
#' @param lfc_threshold Numeric. Log2 fold-change threshold used for coloring.
#'   Default is \code{log2(1.5)}.
#' @param pval_cutoff Numeric. Adjusted p-value cutoff for significance. Default is 0.05.
#' @param min.pct Numeric. min.pvt cutoff for FindMarkers. Default is 0.
#' @param color_palette Named character vector of colors to use for each group.
#'   If NULL, a default palette will be used.
#' @param selected_genes Optional character vector of genes to label on the plots.
#'
#' @details
#' For each row in \code{comparison_table}, the function:
#' \enumerate{
#'   \item Runs \code{FindMarkers()} twice — once for (\code{ident1} vs. \code{ident2})
#'         and once for (\code{ident3} vs. \code{ident4}).
#'   \item Merges results by gene.
#'   \item Assigns each gene to a color group based on the direction and significance
#'         of both comparisons.
#'   \item Plots \eqn{avg_log2FC_1} vs. \eqn{avg_log2FC_2} with color-coded points,
#'         significance thresholds, and labeled selected genes.
#'   \item Saves the resulting plot as a PNG file.
#' }
#'
#' @return Invisibly returns a list of file paths to the saved plots.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- read.csv("2dim_log2fc.csv")
#' user_palette <- c(
#'   "Brown" = "saddlebrown", "Blue" = "dodgerblue", "LightBlue" = "skyblue",
#'   "Green" = "darkgreen", "Purple" = "purple", "Orange" = "orange",
#'   "Pink" = "pink", "LightGreen" = "lightgreen"
#' )
#'
#' double_volcano(object, df, output_dir = "2dim_log2fc",
#'                lfc_threshold = log2(2), pval_cutoff = 0.01,
#'                color_palette = user_palette,
#'                selected_genes = c("CD8A", "PDCD1"))
#' }
double_volcano <- function(
    object,
    comparison_table,
    output_dir = "2dim_log2fc",
    lfc_threshold = log2(1.5),
    pval_cutoff = 0.05,
    min.pct = 0,
    color_palette = NULL,
    selected_genes = NULL
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (is.null(color_palette)) {
    color_palette <- c(
      "Brown" = "darkred", "Blue" = "blue", "LightBlue" = "#27c8e8",
      "Green" = "#008000", "Purple" = "purple", "Orange" = "#FFA500",
      "Pink" = "#ed74d5", "LightGreen" = "#2ff76e"
    )
  }
  
  # Idents(object) <- object$Cohort_tp_mutation
  plot_files <- list()
  merged_list <- list()
  
  for (i in seq_len(nrow(comparison_table))) {
    idents <- unlist(comparison_table[i, 1:4])
    
    markers1 <- Seurat::FindMarkers(object, ident.1 = idents[1], ident.2 = idents[2],
                                    logfc.threshold = 0, min.pct = min.pct)
    markers2 <- Seurat::FindMarkers(object, ident.1 = idents[3], ident.2 = idents[4],
                                    logfc.threshold = 0, min.pct = min.pct)
    
    c1 <- paste0(idents[1], " vs. ", idents[2])
    c2 <- paste0(idents[3], " vs. ", idents[4])
    
    merged <- dplyr::mutate(markers1, gene = rownames(markers1)) %>%
      dplyr::inner_join(
        dplyr::mutate(markers2, gene = rownames(markers2)),
        by = "gene", suffix = c("_1", "_2")
      ) %>%
      dplyr::select(gene, avg_log2FC_1, p_val_adj_1, avg_log2FC_2, p_val_adj_2) %>%
      dplyr::mutate(
        color_group = dplyr::case_when(
          avg_log2FC_1 < -lfc_threshold & avg_log2FC_2 > lfc_threshold &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "Orange",
          avg_log2FC_1 > lfc_threshold & avg_log2FC_2 < -lfc_threshold &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "LightBlue",
          avg_log2FC_1 < -lfc_threshold & avg_log2FC_2 < -lfc_threshold &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "Pink",
          avg_log2FC_1 > lfc_threshold & avg_log2FC_2 > lfc_threshold &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "LightGreen",
          avg_log2FC_1 > lfc_threshold & p_val_adj_1 < pval_cutoff ~ "Blue",
          avg_log2FC_2 > lfc_threshold & p_val_adj_2 < pval_cutoff ~ "Brown",
          avg_log2FC_1 < -lfc_threshold & p_val_adj_1 < pval_cutoff ~ "Green",
          avg_log2FC_2 < -lfc_threshold & p_val_adj_2 < pval_cutoff ~ "Purple",
          TRUE ~ "Other"
        )
      )
    
    # store this merged result with a descriptive name
    merged_list[[paste0(c1, " vs ", c2)]] <- merged
    
    df_gray <- dplyr::filter(merged, color_group == "Other")
    df_color <- dplyr::filter(merged, color_group != "Other")
    
    if (nrow(df_color) == 0) {
      message("No significant genes in row ", i, ". Skipping.")
      next
    }
    
    group_counts <- dplyr::count(df_color, color_group)
    label_positions <- data.frame(
      color_group = c("Brown", "Blue", "Purple", "Green", "LightBlue",
                      "Orange", "LightGreen", "Pink"),
      x = c(0, max(merged$avg_log2FC_1), 0, min(merged$avg_log2FC_1),
            max(merged$avg_log2FC_1), min(merged$avg_log2FC_1),
            max(merged$avg_log2FC_1), min(merged$avg_log2FC_1)),
      y = c(max(merged$avg_log2FC_2), 0, min(merged$avg_log2FC_2), 0,
            min(merged$avg_log2FC_2), max(merged$avg_log2FC_2),
            max(merged$avg_log2FC_2), min(merged$avg_log2FC_2))
    ) %>%
      dplyr::left_join(group_counts, by = "color_group")
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data = df_gray,
                          ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2),
                          color = "gray80", alpha = 0.6, size = 0.8) +
      ggplot2::geom_point(data = df_color,
                          ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2,
                                       color = color_group),
                          alpha = 0.5, size = 0.7) +
      ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                          linetype = "dashed", color = "gray50") +
      ggplot2::geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
                          linetype = "dashed", color = "gray50") +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::labs(x = c1, y = c2, title = "Conserved Log2FC", color = "Group") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(df_color, gene %in% selected_genes),
        ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2, label = gene, color = color_group),
        size = 4, max.overlaps = Inf,
        box.padding = 1, point.padding = 1,
        segment.color = "gray50", segment.size = 0.3,
        arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "closed"),
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = label_positions,
        ggplot2::aes(x = x, y = y, label = n, color = color_group),
        size = 4, fontface = "bold", show.legend = FALSE
      )
    
    outfile <- file.path(output_dir, paste0(gsub(" ", "_", c1), "_", gsub(" ", "_", c2), ".png"))
    ggplot2::ggsave(outfile, p, width = 8, height = 8)
    message("Processed row ", i, ": ", c1, " vs. ", c2)
    print(p)
  }
  
  invisible(merged_list)
}



double_volcano_temp <- function(
    object,
    comparison_table,
    output_dir = "2dim_log2fc",
    lfc_threshold = log2(1.5),
    pval_cutoff = 0.05,
    min.pct = 0,
    color_palette = NULL,
    selected_genes = NULL
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (is.null(color_palette)) {
    color_palette <- c(
      "Brown" = "darkred", "Blue" = "blue", "LightBlue" = "#27c8e8",
      "Green" = "#008000", "Purple" = "purple", "Orange" = "#FFA500",
      "Pink" = "#ed74d5", "LightGreen" = "#2ff76e"
    )
  }
  
  Idents(object) <- object$Cohort_tp_mutation
  plot_files <- list()
  
  for (i in seq_len(nrow(comparison_table))) {
    idents <- unlist(comparison_table[i, 1:4])
    
    markers1 <- Seurat::FindMarkers(object, ident.1 = idents[1], ident.2 = idents[2],
                                    logfc.threshold = 0, min.pct = min.pct)
    markers2 <- Seurat::FindMarkers(object, ident.1 = idents[3], ident.2 = idents[4],
                                    logfc.threshold = 0, min.pct = min.pct)
    
    c1 <- paste0(idents[1], " vs. ", idents[2])
    c2 <- paste0(idents[3], " vs. ", idents[4])
    
    merged <- dplyr::mutate(markers1, gene = rownames(markers1)) %>%
      dplyr::inner_join(
        dplyr::mutate(markers2, gene = rownames(markers2)),
        by = "gene", suffix = c("_1", "_2")
      ) %>%
      dplyr::select(gene, avg_log2FC_1, p_val_adj_1, avg_log2FC_2, p_val_adj_2) %>%
      dplyr::mutate(
        color_group = dplyr::case_when(
          (avg_log2FC_1 < -lfc_threshold | avg_log2FC_2 > lfc_threshold) &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "Orange",
          (avg_log2FC_1 > lfc_threshold | avg_log2FC_2 < -lfc_threshold) &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "LightBlue",
          (avg_log2FC_1 < -lfc_threshold | avg_log2FC_2 < -lfc_threshold) &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "Pink",
          (avg_log2FC_1 > lfc_threshold | avg_log2FC_2 > lfc_threshold) &
            p_val_adj_1 < pval_cutoff & p_val_adj_2 < pval_cutoff ~ "LightGreen",
          avg_log2FC_1 > lfc_threshold & p_val_adj_1 < pval_cutoff ~ "Blue",
          avg_log2FC_2 > lfc_threshold & p_val_adj_2 < pval_cutoff ~ "Brown",
          avg_log2FC_1 < -lfc_threshold & p_val_adj_1 < pval_cutoff ~ "Green",
          avg_log2FC_2 < -lfc_threshold & p_val_adj_2 < pval_cutoff ~ "Purple",
          TRUE ~ "Other"
        )
      )
    
    df_gray <- dplyr::filter(merged, color_group == "Other")
    df_color <- dplyr::filter(merged, color_group != "Other")
    
    if (nrow(df_color) == 0) {
      message("No significant genes in row ", i, ". Skipping.")
      next
    }
    
    group_counts <- dplyr::count(df_color, color_group)
    label_positions <- data.frame(
      color_group = c("Brown", "Blue", "Purple", "Green", "LightBlue",
                      "Orange", "LightGreen", "Pink"),
      x = c(0, max(merged$avg_log2FC_1), 0, min(merged$avg_log2FC_1),
            max(merged$avg_log2FC_1), min(merged$avg_log2FC_1),
            max(merged$avg_log2FC_1), min(merged$avg_log2FC_1)),
      y = c(max(merged$avg_log2FC_2), 0, min(merged$avg_log2FC_2), 0,
            min(merged$avg_log2FC_2), max(merged$avg_log2FC_2),
            max(merged$avg_log2FC_2), min(merged$avg_log2FC_2))
    ) %>%
      dplyr::left_join(group_counts, by = "color_group")
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(data = df_gray,
                          ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2),
                          color = "gray80", alpha = 0.6, size = 0.8) +
      ggplot2::geom_point(data = df_color,
                          ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2,
                                       color = color_group),
                          alpha = 0.5, size = 0.7) +
      ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                          linetype = "dashed", color = "gray50") +
      ggplot2::geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
                          linetype = "dashed", color = "gray50") +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::labs(x = c1, y = c2, title = "Conserved Log2FC", color = "Group") +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(df_color, gene %in% selected_genes),
        ggplot2::aes(x = avg_log2FC_1, y = avg_log2FC_2, label = gene, color = color_group),
        size = 4, max.overlaps = Inf,
        box.padding = 1, point.padding = 1,
        segment.color = "gray50", segment.size = 0.3,
        arrow = grid::arrow(length = grid::unit(0.01, "npc"), type = "closed"),
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = label_positions,
        ggplot2::aes(x = x, y = y, label = n, color = color_group),
        size = 4, fontface = "bold", show.legend = FALSE
      )
    
    outfile <- file.path(output_dir, paste0(gsub(" ", "_", c1), "_", gsub(" ", "_", c2), ".png"))
    ggplot2::ggsave(outfile, p, width = 8, height = 8)
    message("Processed row ", i, ": ", c1, " vs. ", c2)
    plot_files[[i]] <- outfile
  }
  
  invisible(plot_files)
  return(merged)
}

#' Generate Treemaps of TCR Clonotypes for Selected Epitopes
#'
#' This function creates treemaps of TCR clonotypes for each patient and selected epitope,
#' highlighting the distribution of clonotypes. Only patients with a minimum number
#' of cells per epitope are included.
#'
#' @param object A \code{Seurat} object containing TCR metadata in \code{meta.data}.
#' @param top_epitope Character vector of epitopes of interest.
#' @param min_cells Integer. Minimum number of cells for a patient-epitope pair to be included. Default is 20.
#' @param output_dir Character. Directory to save treemap PNGs. Defaults to current working directory.
#' @param color_palette Character vector of colors to use for clonotypes. If NULL, defaults to RColorBrewer "Set3".
#'
#' @return Invisibly returns a list of generated PNG file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' tcr_treemaps(object, top_epitope = c("Epitope1", "Epitope2"), min_cells = 20,
#'                        output_dir = "treemaps/")
#' }
tcr_treemaps <- function(
    object,
    top_epitope,
    min_cells = 20,
    output_dir = ".",
    color_palette = NULL
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  meta_tcr <- object@meta.data
  meta_tcr <- meta_tcr[!is.na(meta_tcr$raw_clonotype_id), ]
  
  # Store valid epitope-patient pairs
  plot_pairs <- data.frame()
  
  for (epitope in top_epitope) {
    temp_meta <- subset(meta_tcr, Specificity == epitope)
    patients <- unique(temp_meta$participant)
    for (p in patients) {
      temp_meta_patient <- subset(temp_meta, participant == p)
      if (nrow(temp_meta_patient) >= min_cells) {
        plot_pairs <- rbind(plot_pairs, 
                            data.frame(epitope = epitope, 
                                       patient = p,
                                       stringsAsFactors = FALSE))
      }
    }
  }
  
  plot_files <- c()
  
  # Process each epitope
  for (epitope in unique(plot_pairs$epitope)) {
    # Generate colors once per epitope
    temp_df <- subset(meta_tcr, Specificity == epitope)
    color_df <- as.data.frame(sort(table(temp_df$raw_clonotype_id), decreasing = TRUE))
    color_df$Freq <- NULL
    colnames(color_df) <- "Clonotype"
    
    if (is.null(color_palette)) {
      base_colors <- RColorBrewer::brewer.pal(12, "Set3")
    } else {
      base_colors <- color_palette
    }
    
    set.seed(match(epitope, unique(plot_pairs$epitope)))
    color_df$Color <- sample(base_colors[(1:nrow(color_df)) %% length(base_colors) + 1])
    
    # Get patients for this epitope
    patients_for_epitope <- plot_pairs$patient[plot_pairs$epitope == epitope]
    
    # Create treemap for each patient
    for (patient in patients_for_epitope) {
      temp_df <- subset(meta_tcr, 
                        participant == patient & Specificity == epitope)
      tree <- left_join(as.data.frame(table(temp_df$raw_clonotype_id)), 
                        color_df,
                        by = c("Var1" = "Clonotype"))
      colnames(tree) <- c("Clonotype", "Count", "Color")
      tree$Clonotype <- paste(tree$Clonotype, paste0("[", tree$Count, "]"))
      
      title_text <- paste(patient, "-", epitope)
      outfile <- file.path(output_dir, paste0(title_text, ".png"))
      
      png(outfile, width = 2000, height = 2000, res = 400, units = "px")
      treemap::treemap(tree,
                       index = "Clonotype",
                       vSize = "Count",
                       vColor = "Color",
                       type = "color",
                       title = title_text,
                       fontsize.labels = 15,
                       border.lwds = 0.5,
                       fontface.labels = 1)
      dev.off()
      
      plot_files <- c(plot_files, outfile)
    }
  }
  
  message("Created ", length(plot_files), " treemap plots in ", output_dir)
  return(invisible(plot_files))
}

#' Generate Patient-Level Expression Heatmap with Z-scores
#'
#' This function creates a heatmap showing gene expression z-scores aggregated
#' at the patient level. Expression values are first averaged per patient within
#' each cohort/timepoint, then z-score normalized across all samples.
#'
#' @param object A Seurat object containing expression data and metadata
#' @param features Character vector of gene names to include in the heatmap
#' @param groups Character vector specifying the order of cohort/timepoint groups.
#'   This defines both which groups to include and their display order.
#' @param idents Character string specifying the metadata column containing
#'   cohort/timepoint information. If NULL, defaults to "cohort_tp".
#' @param pseudobulk_column Character string specifying the metadata column containing
#'   pseudobulk. If NULL, defaults to "Patient".
#' @param colors Named vector of colors for each cohort/timepoint group. Names must
#'   match the values in the idents If NULL, colors will be automatically
#'   generated from a standard palette.
#' @param cluster_rows Logical indicating whether to cluster rows (genes).
#'   Default is FALSE.
#' @param cluster_columns Logical indicating whether to cluster columns (samples).
#'   Default is FALSE.
#' @param show_column_names Logical indicating whether to show sample names.
#'   Default is FALSE.
#' @param heatmap_title Character string for the heatmap legend title.
#'   Default is "Z-score".
#' @param color_range Numeric vector of length 3 specifying the values for
#'   (min, center, max) of the color scale. Default is c(-2, 0, 2).
#' @param color_palette Character vector of length 3 specifying colors for
#'   (min, center, max). Default is c("#4575b4", "white", "#d73027").
#' @param annotation_height Unit object specifying the height of the annotation bar.
#'   Default is unit(4, "mm").
#'
#' @return A ComplexHeatmap object
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts normalized expression data for specified genes
#'   \item Averages expression per gene, per patient, per cohort/timepoint
#'   \item Computes z-scores across all samples for each gene (row-wise)
#'   \item Orders samples according to the specified group order
#'   \item Creates a heatmap with cohort annotations
#' }
#'
#' Z-scores are calculated as: \deqn{z = (x - mean(x)) / sd(x)}
#' where x represents expression values across all samples for a given gene.
#' NA values (from genes with constant expression) are replaced with 0.
#'
#' If colors are not provided, the function automatically generates colors using
#' a standard palette based on the unique values in the cohort column.
#'
#' @examples
#' \dontrun{
#' # Basic usage with auto-generated colors
#' heatmap_pseudobulk(
#'   object = object,
#'   features = c("CD3D", "CD8A", "IFNG"),
#'   groups = c("Chronic Pre-DAA", "ACTG Pre-DAA", "SR Pre-resolution")
#' )
#'
#' # With custom colors
#' my_colors <- c(
#'   "Chronic Pre-DAA" = "#F8A19FB2",
#'   "ACTG Pre-DAA" = "#2ED9FFB2",
#'   "SR Pre-resolution" = "#FEAF16B2"
#' )
#'
#' heatmap_pseudobulk(
#'   object = object,
#'   features = marker_genes,
#'   groups = c("Chronic Pre-DAA", "ACTG Pre-DAA", "SR Pre-resolution"),
#'   colors = my_colors,
#'   cluster_rows = TRUE
#' )
#' 
#' # With custom metadata columns
#' heatmap_pseudobulk(
#'   object = object,
#'   features = marker_genes,
#'   groups = c("Group1", "Group2", "Group3"),
#'   idents = "treatment_group",
#'   pseudobulk_column = "sample_id"
#' )
#' }
#'
#' @export
heatmap_pseudobulk <- function(object,
                                            features,
                                            groups,
                                            idents = NULL,
                                            pseudobulk_column = NULL,
                                            colors = NULL,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            show_column_names = FALSE,
                                            heatmap_title = "Z-score",
                                            color_range = c(-2, 0, 2),
                                            color_palette = c("#4575b4", "white", "#d73027"),
                                            annotation_height = unit(4, "mm")) {
  
  # Set default column names if NULL
  if (is.null(idents)) {
    idents <- "cohort_tp"
  }
  
  if (is.null(pseudobulk_column)) {
    pseudobulk_column <- "Patient"
  }
  
  cells_to_keep <- colnames(object)[object[[idents]] %>% unlist %in% groups]
  object <- subset(object, cells = cells_to_keep)
  
  features <- intersect(features, rownames(object))
  
  # --- Extract expression data ---
  expression_matrix <- as.data.frame(GetAssayData(object, layer = "data"))
  expression_matrix <- expression_matrix[features, ]
  
  # --- Metadata for cohort and patient ---
  cell_metadata <- data.frame(
    cell = colnames(expression_matrix),
    cohort_tp = object[[idents, drop = TRUE]],
    patient_id = object[[pseudobulk_column, drop = TRUE]]
  )
  
  # --- Aggregate expression by patient ---
  grouped_expression <- expression_matrix %>%
    tibble::rownames_to_column(var = "gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(cell_metadata, by = "cell") %>%
    group_by(gene, cohort_tp, patient_id) %>%
    summarize(mean_expression = mean(expression), .groups = "drop") %>%
    mutate(sample_id = paste(cohort_tp, patient_id, sep = " ")) %>%
    select(gene, sample_id, mean_expression) %>%
    pivot_wider(
      names_from = sample_id,
      values_from = mean_expression,
      values_fill = list(mean_expression = 0)
    )
  
  # --- Convert to matrix ---
  grouped_expression_matrix <- as.matrix(grouped_expression[, -c(1)])
  rownames(grouped_expression_matrix) <- grouped_expression$gene
  grouped_expression_matrix <- grouped_expression_matrix[features, ]
  
  # --- Compute z-scores per gene (row) ---
  zscore_matrix <- t(scale(t(grouped_expression_matrix)))
  zscore_matrix[is.na(zscore_matrix)] <- 0
  
  # --- Prepare cohort grouping with custom order ---
  cohort_labels <- sapply(strsplit(colnames(zscore_matrix), " "), `[`, 1)
  
  # Convert to factor with custom levels
  cohort_labels <- factor(cohort_labels, levels = groups)
  
  # Reorder matrix columns by cohort_labels
  ord <- order(cohort_labels)
  zscore_matrix <- zscore_matrix[, ord]
  cohort_labels <- cohort_labels[ord]
  unique_cohorts <- unique(as.character(cohort_labels))
  # Generate colors if not provided
  if (is.null(colors)) {
    n_colors <- length(unique_cohorts)
    if (n_colors <= 8) {
      palette_colors <- RColorBrewer::brewer.pal(max(3, n_colors), "Set2")
    } else if (n_colors <= 12) {
      palette_colors <- RColorBrewer::brewer.pal(n_colors, "Set3")
    } else {
      palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_colors)
    }
    colors <- scales::alpha(palette_colors[1:n_colors], 0.8)
  }else if(length(colors) < length(unique_cohorts)){
    colors <- rep(colors, length.out = length(unique_cohorts))
  }
  colors <- colors[1: length(unique_cohorts)]
  names(colors) <- unique_cohorts
  
  
  # --- Column annotation ---
  ha <- HeatmapAnnotation(
    Cohort = cohort_labels,
    col = list(Cohort = colors),
    show_legend = TRUE,
    annotation_name_side = "left",
    simple_anno_size = annotation_height
  )
  
  # --- Draw ComplexHeatmap ---
  hm <- Heatmap(
    zscore_matrix,
    name = heatmap_title,
    top_annotation = ha,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_column_names = show_column_names,
    column_split = cohort_labels,
    col = circlize::colorRamp2(color_range, color_palette),
    heatmap_legend_param = list(title = heatmap_title)
  )
  
  return(hm)
}


#' Generate Cell-Type-Level Expression Heatmap
#'
#' This function creates a heatmap showing gene expression averaged across cell types
#' or groups. Expression values are aggregated by taking the mean across all cells
#' within each group, then optionally row-scaled (z-score normalized) for visualization.
#' It allows for customization of gene label appearance.
#'
#' @param object A Seurat object containing expression data and metadata
#' @param features Character vector of gene names to include in the heatmap
#' @param groups Character vector specifying the order of cell type/group labels.
#' This defines both which groups to include and their display order.
#' @param idents Character string specifying the metadata column containing
#' cell type or group information. If NULL, defaults to "cohort_tp".
#' @param bold_genes Character vector of gene names to render in **bold** text.
#' Uses \code{bquote(bold(...))} for formatting. Default is NULL (no bolding).
#' @param row_fontsize Numeric value specifying the font size for row labels (gene names).
#' Default is 10.
#' @param normalize Logical indicating whether to normalize the data before plotting.
#' Default is TRUE.
#' @param cluster_rows Logical indicating whether to cluster rows (genes).
#' Default is FALSE.
#' @param cluster_cols Logical indicating whether to cluster columns (groups).
#' Default is FALSE.
#' @param scale Character string indicating scaling method. Options are "row" (default),
#' "column", or "none".
#' @param title Optional character string for plot title. Default is NA (no title).
#' @param angle_col Character string or numeric specifying the angle for column labels.
#' Default is "315".
#' @param file_name name of png file.
#' @param ... Additional arguments passed to pheatmap::pheatmap()
#'
#' @return A pheatmap object
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#' \item Optionally normalizes the Seurat object
#' \item Extracts expression data for specified genes
#' \item Averages expression across all cells within each group
#' \item Orders groups according to the specified order
#' \item Creates a heatmap with row scaling (z-score by default), applying bolding
#' to specified gene names via \code{labels_row} and setting global font size via
#' \code{fontsize_row}.
#' }
#'
#' Row scaling computes z-scores as: $z = (x - mean(x)) / sd(x)$
#' where x represents expression values across all groups for a given gene.
#'
#' @examples
#' \dontrun{
#' # Basic usage with bolded genes and larger font size
#' heatmap_meta(
#' object = object,
#' features = c("CD3D", "CD8A", "IFNG", "FOXP3"),
#' groups = c("CD4 T", "CD8 T", "NK"),
#' bold_genes = c("CD8A", "FOXP3"),
#' row_fontsize = 12
#' )
#'
#' # With custom group column and no normalization
#' heatmap_meta(
#' object = object,
#' features = marker_genes,
#' groups = c("Chronic Pre-DAA", "ACTG Pre-DAA", "SR Pre-resolution"),
#' idents = "condition",
#' normalize = FALSE
#' )
#' }
#'
#' @export
heatmap_meta <- function(object,
                         features,
                         groups,
                         idents = NULL,
                         bold_genes = NULL, 
                         row_fontsize = 10,
                         normalize = TRUE,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         scale = "row",
                         title = NA,
                         angle_col = "315",
                         file_name = NULL,
                         ...) {
  
  # Set default column name if NULL
  if (is.null(idents)) {
    idents <- "cohort_tp"
  }
  
  # Normalize data if requested (using Seurat's NormalizeData)
  if (normalize) {
    object <- Seurat::NormalizeData(object)
  }
  
  features <- intersect(features, rownames(object))
  
  # --- Extract expression data ---
  expression_matrix <- as.data.frame(Seurat::GetAssayData(object, layer = "data"))
  expression_matrix <- expression_matrix[features, ]
  
  # --- Add cell type annotations ---
  cell_annotations <- object[[idents, drop = TRUE]]
  
  # --- Summarize the expression data by cell type (mean) ---
  grouped_expression <- expression_matrix %>%
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    dplyr::left_join(
      data.frame(cell = colnames(expression_matrix), cell_type = cell_annotations),
      by = "cell"
    ) %>%
    dplyr::group_by(gene, cell_type) %>%
    dplyr::summarize(mean_expression = mean(expression), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = mean_expression)
  
  # --- Convert to matrix and order rows/columns ---
  grouped_expression_matrix <- as.matrix(grouped_expression[, -1])
  rownames(grouped_expression_matrix) <- grouped_expression$gene
  # Ensure the order matches the input 'features'
  grouped_expression_matrix <- grouped_expression_matrix[features, ]
  # Ensure the column order matches the input 'groups'
  grouped_expression_matrix <- grouped_expression_matrix[, groups]
  
  # --- **Logic to create custom row labels (gene names)** ---
  if (!is.null(bold_genes) && any(bold_genes %in% rownames(grouped_expression_matrix))) {
    
    # Create an expression vector. The 'bquote' function is used to create an 
    # expression object that pheatmap can interpret for font styles (bold).
    labels_expression <- lapply(rownames(grouped_expression_matrix), function(gene) {
      if (gene %in% bold_genes) {
        # Format the gene name as bold using bquote(bold(...))
        return(bquote(bold(.(as.character(gene)))))
      } else {
        # Keep other gene names as standard text
        return(bquote(.(as.character(gene))))
      }
    })
    
    # Convert the list of expressions into a single expression vector
    labels_expression <- as.expression(labels_expression)
  } else {
    # If no bold genes are specified, use the original row names
    labels_expression <- rownames(grouped_expression_matrix)
  }
  
  # --- Create heatmap ---
  p <- pheatmap::pheatmap(
    grouped_expression_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale = scale,
    main = title,
    angle_col = angle_col,
    fontsize_row = row_fontsize,
    labels_row = labels_expression, 
    ...
  )
  
  if(!is.null(file_name)){
    width <- length(unique(groups)) * 0.5 + 5
    height <- length(unique(features)) * 0.35 + 6
    png(filename = file_name, 
        width = width, 
        height = height,
        units = "cm",
        res = 200)
    print(p)
    dev.off()
  }
  
  return(p)
}



#' Create Venn Diagram from Differential Expression Comparisons
#'
#' This function performs multiple pairwise differential expression comparisons
#' between groups and visualizes the overlap of significantly differentially
#' expressed genes using an area-proportional Venn/Euler diagram.
#'
#' @param object A Seurat object containing expression data and metadata
#' @param comparisons Named list where names are comparison labels
#'   (e.g., "Group A vs Group B") and values are character vectors of length 2
#'   specifying which groups to compare (e.g., c("Group A", "Group B")).
#'   Must contain between 3 and 6 comparisons.
#' @param idents Character string specifying the metadata column containing
#'   group identities. This column will be set as the active identity.
#' @param direction Character string specifying which genes to include in the Venn.
#'   Options are "up" (upregulated in ident.1), "down" (downregulated in ident.1),
#'   or "both". Default is "up".
#' @param log2fc_threshold Numeric value specifying the minimum absolute log2
#'   fold change threshold for significance. Default is 0.58.
#' @param min_pct Numeric value specifying the minimum fraction of cells expressing
#'   the gene in either group. Default is 0.3.
#' @param p_adj_cutoff Numeric value for adjusted p-value cutoff. Default is 0.05.
#' @param fills Character vector of colors for the Venn diagram circles. If NULL,
#'   a default palette will be used. Length should match the number of comparisons.
#' @param alpha Numeric value between 0 and 1 for transparency of circles.
#'   Default is 0.6.
#' @param main Character string for the main title. If NULL, a default title
#'   will be generated. Default is NULL.
#' @param subtitle Character string for subtitle. If NULL, parameters will be
#'   shown. Default is NULL.
#' @param label_cex Numeric value for label text size. Default is 1.2.
#' @param quantity_cex Numeric value for quantity text size. Default is 1.0.
#' @param file_name Output png file name. Default is eulerr_plot.
#' @param ... Additional arguments passed to eulerr::euler()
#'
#' @return A list containing:
#' \describe{
#'   \item{markers}{A named list of data frames, one per comparison, containing
#'     full differential expression results from FindMarkers}
#'   \item{genes}{A named list of character vectors containing the filtered gene
#'     names that passed the thresholds for each comparison}
#'   \item{euler_fit}{The fitted euler object}
#'   \item{plot}{The plot is drawn directly; returns invisible}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Sets the specified metadata column as active identity
#'   \item Performs pairwise differential expression using Seurat::FindMarkers
#'   \item Filters genes based on p_adj, log2FC, and direction criteria
#'   \item Creates an area-proportional Venn/Euler diagram showing overlaps
#' }
#'
#' For "up" direction, genes are filtered where avg_log2FC > log2fc_threshold.
#' For "down" direction, genes are filtered where avg_log2FC < -log2fc_threshold.
#' For "both" direction, genes are filtered where |avg_log2FC| > log2fc_threshold.
#'
#' @examples
#' \dontrun{
#' # Custom comparisons with named list
#' comparisons <- list(
#'   "Chronic vs ACTG" = c("Chronic Post-DAA Conserved", "ACTG Post-DAA Conserved"),
#'   "Chronic vs SR" = c("Chronic Post-DAA Conserved", "SR ≥24W-post-resolution unknown"),
#'   "ACTG vs SR" = c("ACTG Post-DAA Conserved", "SR ≥24W-post-resolution unknown")
#' )
#'
#' result <- venn_plots(
#'   object = object,
#'   comparisons = comparisons,
#'   idents = "Cohort_tp_mutation",
#'   log2fc_threshold = 0.58,
#'   min_pct = 0.3,
#'   main = "Upregulated Genes"
#' )
#'
#' # Look at downregulated genes
#' comparisons <- list(
#'   "A vs B" = c("Group1", "Group2"),
#'   "A vs C" = c("Group1", "Group3"),
#'   "B vs C" = c("Group2", "Group3")
#' )
#'
#' venn_plots(
#'   object = object,
#'   comparisons = comparisons,
#'   idents = "treatment",
#'   direction = "down",
#'   log2fc_threshold = 0.5
#' )
#'
#' # Access results
#' result$markers  # Full marker data frames
#' result$genes    # Filtered gene lists
#' }
#'
#' @export
venn_plots <- function(object,
                              comparisons,
                              idents,
                              direction = "up",
                              log2fc_threshold = 0.58,
                              min_pct = 0.3,
                              p_adj_cutoff = 0.05,
                              fills = NULL,
                              alpha = 0.6,
                              main = NULL,
                              subtitle = NULL,
                              label_cex = 1.2,
                              quantity_cex = 1.0,
                              file_name = "eulerr_plot",
                              ...) {
  
  # Input validation
  if (!is.list(comparisons) || is.null(names(comparisons))) {
    stop("comparisons must be a named list")
  }
  
  if (length(comparisons) < 3 || length(comparisons) > 6) {
    stop("Number of comparisons must be between 3 and 6")
  }
  
  if (!idents %in% colnames(object@meta.data)) {
    stop("Group column '", idents, "' not found in metadata")
  }
  
  if (!direction %in% c("up", "down", "both")) {
    stop("Direction must be 'up', 'down', or 'both'")
  }
  
  # Validate each comparison
  for (comp_name in names(comparisons)) {
    comp <- comparisons[[comp_name]]
    if (length(comp) != 2) {
      stop("Each comparison must contain exactly 2 groups. Check: ", comp_name)
    }
  }
  
  # Set active identity
  Idents(object) <- object[[idents]] %>% unlist() %>% as.character()
  
  # Perform differential expression for each comparison
  markers_list <- list()
  genes_list <- list()
  
  message("Performing differential expression comparisons...")
  for (comp_name in names(comparisons)) {
    comp <- comparisons[[comp_name]]
    
    message("  ", comp_name, ": ", comp[1], " vs ", comp[2])
    
    markers <- FindMarkers(
      object,
      ident.1 = comp[1],
      ident.2 = comp[2],
      logfc.threshold = 0,
      min.pct = min_pct
    )
    
    markers$gene <- rownames(markers)
    markers_list[[comp_name]] <- markers
    
    # Filter genes based on direction
    if (direction == "up") {
      filtered_genes <- rownames(
        markers %>% filter(p_val_adj < p_adj_cutoff, avg_log2FC > log2fc_threshold)
      )
    } else if (direction == "down") {
      filtered_genes <- rownames(
        markers %>% filter(p_val_adj < p_adj_cutoff, avg_log2FC < -log2fc_threshold)
      )
    } else {  # both
      filtered_genes <- rownames(
        markers %>% filter(p_val_adj < p_adj_cutoff, abs(avg_log2FC) > log2fc_threshold)
      )
    }
    
    genes_list[[comp_name]] <- filtered_genes
    message("    Found ", length(filtered_genes), " genes")
  }
  
  # Default colors if not provided
  if (is.null(fills)) {
    n_comp <- length(comparisons)
    if (n_comp <= 3) {
      fills <- c("#66c2a5", "#fc8d62", "#8da0cb")
    } else if (n_comp <= 6) {
      fills <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f")
    } else {
      fills <- rainbow(n_comp)
    }
  }
  
  # Compute Euler fit
  fit <- euler(genes_list, ...)
  
  # Generate main title if not provided
  if (is.null(main)) {
    direction_text <- switch(direction,
                             "up" = "Upregulated",
                             "down" = "Downregulated",
                             "both" = "Differentially Expressed")
    main <- paste(direction_text, "Genes")
  }
  
  # Generate subtitle if not provided
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "min.pct = ", min_pct, 
      "  |log2FC| > ", log2fc_threshold,
      "  p.adj < ", p_adj_cutoff
    )
  }
  
  # Create plot
  png(file_name, width = 2000, height = 1600, res = 300)
  print(plot(
    fit,
    fills = fills[1:length(genes_list)],
    alpha = alpha,
    labels = list(font = 1, cex = label_cex),
    quantities = list(cex = quantity_cex),
    edges = TRUE,
    main = main
  ))
  # Add subtitle
  grid.text(
    subtitle,
    y = unit(0.93, "npc"),
    gp = gpar(cex = 0.8)
  )
  dev.off()
  
  # Return results invisibly
  invisible(list(
    markers = markers_list,
    genes = genes_list,
    euler_fit = fit
  ))
}


#' Plot Top N Usage Columns From cNMF Gene Spectra
#'
#' This function reads cNMF gene spectra, normalizes them, extracts gene-level
#' statistics (min/median/max) from a Seurat object, and produces barplots for
#' the top-N contributing genes per program. It also returns a top-gene table.
#'
#' @param dir_path Character. Directory containing cNMF output files.
#' @param object Seurat or SingleCellExperiment object containing normalized data.
#' @param output_dir Directory where PNG plots should be saved.
#' @param top_n Number of top genes to extract and plot (default = 30).
#' @param show_metrics Logical. If TRUE, add min/median/max labels to bars.
#' @param highlight_genes A gene list to be bold and red on barplots. If NULL, output without hightlight.
#'
#' @return A list containing:
#'   \item{top_df}{Data frame of top-N genes per program}
#'   \item{file_used}{Path to the CNMF spectrum file used}
#'
#' @export
cnmf_bar_plots <- function(
    dir_path,
    object,
    output_dir,
    top_n = 30,
    show_metrics = FALSE,
    highlight_genes = NULL
) {
  
  # --- Setup and Data Loading (Unchanged) ---
  
  if (!dir.exists(dir_path))
    stop("dir_path does not exist.")
  
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  
  file_name <- list.files(
    dir_path,
    pattern = "cnmf_run\\.gene_spectra_tpm.*\\.dt_0_20\\.txt$",
    full.names = TRUE
  )
  if (length(file_name) == 0)
    stop("No CNMF spectrum files found.")
  
  file_name <- file_name[1]
  gene_by_pro <- read.table(file_name, fill = TRUE, header = TRUE, row.names = 1)
  
  # --- Extract top-N genes (table output) ---
  
  # The top_df will store the raw gene names
  top_matrix <- matrix(nrow = nrow(gene_by_pro), ncol = top_n)
  rownames(top_matrix) <- rownames(gene_by_pro)
  
  for (i in seq_len(nrow(gene_by_pro))) {
    row_data <- as.numeric(gene_by_pro[i, ])
    top_idx <- order(row_data, decreasing = TRUE)[1:top_n]
    top_matrix[i, ] <- colnames(gene_by_pro)[top_idx]
  }
  
  top_df <- as.data.frame(t(top_matrix))
  colnames(top_df) <- paste0("Usage_", colnames(top_df))
  
  # --- Prepare output data frame for highlighting ---
  if (!is.null(highlight_genes)) {
    # Function to apply highlighting format (Markdown/HTML)
    format_gene <- function(gene) {
      if (gene %in% highlight_genes) {
        # Use HTML/Markdown for bold and red text
        return(paste0("<span style='color:red;'>**", gene, "**</span>"))
      } else {
        return(gene)
      }
    }
  }
  
  # --- Normalize data and extract metrics ---
  gene_by_pro <- t(apply(gene_by_pro, 1, function(x) (x - min(x)) / (max(x) - min(x))))
  
  norm_data <- SeuratObject::GetAssayData(object, layer = "data") |> as.matrix()
  norm_data <- norm_data[colnames(gene_by_pro), ]
  
  gene_medians <- matrixStats::rowMedians(norm_data)
  gene_maxs    <- matrixStats::rowMaxs(norm_data)
  
  norm_data[norm_data == 0] <- NA
  gene_mins <- matrixStats::rowMins(norm_data, na.rm = TRUE)
  
  # Dynamic plot size (Unchanged)
  width_px  <- 1000
  height_px <- 1000 + ((top_n - 30) * 20)
  
  # --- Helper Plot Function (Modified) ---
  plot_top_columns <- function(row_data, row_name, gene_mins, gene_medians, gene_maxs, highlight_genes) {
    
    top_idx    <- order(row_data, decreasing = TRUE)[1:top_n] |> rev()
    top_vals   <- row_data[top_idx] |> rev()
    top_genes  <- names(top_vals)
    
    # Apply highlighting format for PLOT LABELS
    formatted_genes <- sapply(top_genes, function(gene) {
      if (!is.null(highlight_genes) && gene %in% highlight_genes) {
        # Use Markdown/HTML syntax for ggtext
        return(paste0("<span style='color:red;'>**", gene, "**</span>"))
      } else {
        return(gene)
      }
    })
    
    # Create factors from formatted genes
    formatted_genes_factor <- factor(formatted_genes, levels = rev(formatted_genes))
    
    if (show_metrics) {
      plot_data <- data.frame(
        Column     = formatted_genes_factor,
        Value      = top_vals,
        MinExpr    = round(gene_mins[top_genes], 4),
        MedianExpr = round(gene_medians[top_genes], 4),
        MaxExpr    = round(gene_maxs[top_genes], 4)
      )
    } else {
      plot_data <- data.frame(
        Column = formatted_genes_factor,
        Value  = top_vals
      )
    }
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Column, y = Value, fill = Value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(
        title = paste("Usage", row_name),
        x = "",
        y = "Values"
      ) +
      ggplot2::theme_classic() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = "#68bdde", high = "#de6868") +
      # Use ggtext::element_markdown() to interpret the HTML/Markdown in the labels
      ggplot2::theme(
        axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
      )
    
    if (show_metrics) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = paste0(MinExpr, "  ", MedianExpr, "  ", MaxExpr)),
          hjust = 1,
          size = 2
        )
    }
    
    outfile <- file.path(output_dir, paste0("Usage_", row_name, ".png"))
    grDevices::png(outfile, width = width_px, height = height_px, res = 200)
    print(p)
    grDevices::dev.off()
  }
  
  # --- Loop and plot ---
  for (i in seq_len(nrow(gene_by_pro))) {
    row_name <- rownames(gene_by_pro)[i]
    row_data <- gene_by_pro[i, ]
    # Pass the highlight_genes argument to the helper function
    plot_top_columns(row_data, row_name, gene_mins, gene_medians, gene_maxs, highlight_genes)
  }
  
  message("Finished generating plots in: ", output_dir)
  
  return(list(
    # Return the data frame with HTML/Markdown formatting
    top_df = top_df, 
    file_used = file_name
  ))
}


#' Plot UMAP Usage Maps From cNMF Usage Matrix
#'
#' This function reads the cNMF "usages" consensus matrix, normalizes usage
#' values per cell, merges them into the Seurat object's metadata, and produces
#' UMAP feature plots for all Usage_* programs.
#'
#' @param dir_path Directory containing the cNMF usages file.
#' @param object Seurat object with UMAP reduction.
#' @param output_dir Directory where UMAP PNG files should be saved.
#' @param pt_size Point size for FeaturePlot (default = 0.5).
#'
#' @return The updated Seurat object (with usage values added to metadata).
#' @export
cnmf_umaps <- function(
    dir_path,
    object,
    output_dir,
    pt_size = 0.5
) {
  
  # Validate Inputs
  if (!dir.exists(dir_path))
    stop("dir_path does not exist.")
  
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)
  
  # Locate Usage Matrix
  file_name <- list.files(
    dir_path,
    pattern = "cnmf_run\\.usages.*\\.dt_0_20\\.consensus\\.txt$",
    full.names = TRUE
  )
  if (length(file_name) == 0)
    stop("No cNMF usage consensus file found in dir_path.")
  
  file_name <- file_name[1]
  
  # Read + Normalize Usage Matrix
  cell_by_pro <- read.table(file_name, fill = TRUE, header = TRUE, row.names = 1)
  
  # Normalize per-cell (row sums = 1)
  cell_by_pro <- t(apply(cell_by_pro, 1, function(x) x / sum(x)))
  cell_by_pro <- as.data.frame(cell_by_pro)
  
  # Clean column names
  colnames(cell_by_pro) <- gsub("^X", "Usage", colnames(cell_by_pro))
  cell_by_pro$cell_id <- rownames(cell_by_pro)
  
  # Add Usage to Seurat Metadata
  meta <- object@meta.data
  meta$cell_id <- rownames(meta)
  
  meta <- dplyr::inner_join(meta, cell_by_pro, by = "cell_id")
  rownames(meta) <- meta$cell_id
  
  object@meta.data <- meta
  
  # Identify Usage Features
  usage_features <- grep("^Usage", colnames(object@meta.data), value = TRUE)
  
  # Generate UMAP Plots
  for (usage in usage_features) {
    outfile <- file.path(output_dir, paste0("umap_", usage, ".png"))
    grDevices::png(outfile, width = 1000, height = 1000, res = 200)
    
    FeaturePlot_scCustom(
      object,
      features = usage,
      pt.size = pt_size,
      order = TRUE
    ) |> print()
    
    grDevices::dev.off()
  }
  
  message("UMAP usage plots saved to: ", output_dir)
  
  return(object)
}


#' Search Gene Ranks Across cNMF Programs
#'
#' @param dir_path Path to the directory containing cNMF output files
#' @param features Character vector of gene names to search for
#' @return A data frame with genes as rows and programs as columns, showing ranks
#' @export
cnmf_gene_ranks <- function(dir_path, features) {
  
  # --- Validate inputs ---
  if (!dir.exists(dir_path)) {
    stop("dir_path does not exist.")
  }
  
  if (!is.character(features) || length(features) == 0) {
    stop("genes must be a non-empty character vector.")
  }
  
  # --- Load cNMF gene spectra file ---
  file_name <- list.files(
    dir_path,
    pattern = "cnmf_run\\.gene_spectra_tpm.*\\.dt_0_20\\.txt$",
    full.names = TRUE
  )
  
  if (length(file_name) == 0) {
    stop("No CNMF spectrum files found.")
  }
  
  file_name <- file_name[1]
  gene_by_pro <- read.table(file_name, fill = TRUE, header = TRUE, row.names = 1)
  
  # --- Check which genes are present ---
  genes_present <- features[features %in% colnames(gene_by_pro)]
  genes_missing <- features[!features %in% colnames(gene_by_pro)]
  
  if (length(genes_missing) > 0) {
    warning("The following genes were not found in cNMF results: ", 
            paste(genes_missing, collapse = ", "))
  }
  
  if (length(genes_present) == 0) {
    stop("None of the specified genes were found in the cNMF results.")
  }
  
  # --- Calculate ranks for each program ---
  rank_matrix <- matrix(
    nrow = length(genes_present),
    ncol = nrow(gene_by_pro),
    dimnames = list(genes_present, rownames(gene_by_pro))
  )
  
  for (i in seq_len(nrow(gene_by_pro))) {
    program_name <- rownames(gene_by_pro)[i]
    row_data <- as.numeric(gene_by_pro[i, ])
    names(row_data) <- colnames(gene_by_pro)
    
    # Rank genes (1 = highest weight)
    ranks <- rank(-row_data, ties.method = "min")
    
    # Extract ranks for genes of interest
    rank_matrix[, program_name] <- ranks[genes_present]
  }
  
  # --- Format output ---
  result <- as.data.frame(rank_matrix)
  colnames(result) <- paste0("Usage_", colnames(result))
  
  cat("Total genes in analysis:", ncol(gene_by_pro), "\n")
  cat("Total programs in analysis:", nrow(gene_by_pro), "\n")
  cat("======================\n")
  print(result)
  
  return(result)
}



#' Add top cNMF programs to Seurat metadata with usage cutoff
#'
#' This function loads a cNMF consensus usage file, normalizes usage values,
#' identifies the top 4 programs per cell (Primary–Quaternary) above a cutoff,
#' and adds them to the metadata of a Seurat object. Programs with usage
#' values below or equal to the cutoff are ignored and filled with NA.
#'
#' @param obj A Seurat object.
#' @param dir_path Path to the directory containing the cNMF usage file.
#' @param cutoff Numeric. Minimum normalized usage required for a program
#'        to be labeled. Default is 0.
#'
#' @return A Seurat object with these metadata columns added:
#'   \itemize{
#'     \item Primary
#'     \item Secondary
#'     \item Tietiary
#'     \item Quaternary
#'   }
#'   and a row-normalized cNMF cells x programs matrix (sum in each row = 1).
#'
#' @examples
#' \dontrun{
#' obj <- cnmf_top_programs(obj, "path/to/cnmf/output/", cutoff = 0.03)
#' }
#'
#' @export
cnmf_top_programs <- function(obj,
                              dir_path,
                              cutoff = 0) {
  # --- Remove old assignments ---
  cols_to_remove <- c("Primary", "Secondary", "Tietiary", "Quaternary")
  existing <- intersect(cols_to_remove, colnames(obj@meta.data))
  if (length(existing) > 0) {
    obj@meta.data <- obj@meta.data[, !colnames(obj@meta.data) %in% existing, drop = FALSE]
  }
  
  # --- Locate usage file ---
  file_name <- list.files(
    dir_path,
    pattern = "cnmf_run\\.usages.*consensus\\.txt$",
    full.names = TRUE
  )
  
  if (length(file_name) == 0) {
    stop("No cNMF usage file found in dir_path.")
  }
  
  # Use first match
  file_name <- file_name[1]
  
  # --- Load and normalize ---
  cell_by_pro <- read.table(file_name, fill = TRUE, header = TRUE, row.names = 1)
  
  # Row-normalize each cell's usages
  cell_by_pro <- t(apply(cell_by_pro, 1, function(x) x / sum(x)))
  cell_by_pro <- as.data.frame(cell_by_pro)
  
  colnames(cell_by_pro) <- gsub("^X", "Usage", colnames(cell_by_pro))
  cell_by_pro$cell_id <- rownames(cell_by_pro)
  
  # --- Reshape + apply cutoff ---
  df <- cell_by_pro %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(setdiff(colnames(cell_by_pro), "cell_id")),
      names_to = "usage",
      values_to = "Expression"
    ) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::arrange(desc(Expression), .by_group = TRUE) %>%
    dplyr::mutate(Expression = ifelse(Expression > cutoff, Expression, NA_real_),
                  usage = ifelse(is.na(Expression), NA_character_, usage)) %>%
    dplyr::summarise(
      Primary    = usage[1],
      Secondary  = usage[2],
      Tietiary   = usage[3],
      Quaternary = usage[4],
      .groups = "drop"
    )
  
  # --- Merge with metadata ---
  meta <- obj@meta.data
  meta <- dplyr::left_join(meta, df, by = "cell_id")
  rownames(meta) <- meta$cell_id
  obj@meta.data <- meta
  
  # --- Factor ordering ---
  obj$Primary    <- factor(obj$Primary,    levels = stringr::str_sort(unique(obj$Primary),    numeric = TRUE))
  obj$Secondary  <- factor(obj$Secondary,  levels = stringr::str_sort(unique(obj$Secondary),  numeric = TRUE))
  obj$Tietiary   <- factor(obj$Tietiary,   levels = stringr::str_sort(unique(obj$Tietiary),   numeric = TRUE))
  obj$Quaternary <- factor(obj$Quaternary, levels = stringr::str_sort(unique(obj$Quaternary), numeric = TRUE))
  
  cell_by_pro$cell_id <- NULL
  return(list(object = obj, cell_by_pro = cell_by_pro))
}


#' Compare Overlap Between Two cNMF Program Sets
#'
#' This function automatically finds cNMF *gene_spectra_tpm* files inside two
#' directories, extracts the top-N genes from each program, optionally filters
#' to variable features, and plots an overlap heatmap.
#'
#' @param dir1 Directory containing the first CNMF run.
#' @param dir2 Directory containing the second CNMF run.
#' @param object Seurat object (required if vf_only = TRUE).
#' @param top_n Number of top genes to use per program.
#' @param label1 Label prefix for heatmap rows.
#' @param label2 Label prefix for heatmap columns.
#' @param heatmap_title Title to display on the heatmap.
#' @param vf_only Logical; if TRUE, keep only genes in VariableFeatures(object).
#'
#' @return A list: overlap_matrix, top_df1, top_df2
#' @export
compare_cnmf_programs <- function(
    dir1,
    dir2,
    object,
    top_n = 50,
    label1 = "A",
    label2 = "B",
    heatmap_title = "Usage Similarity",
    vf_only = TRUE
) {
  
  # --- Internal helper: find the CNMF spectra file ---
  find_cnmf_file <- function(dir_path) {
    file <- list.files(
      dir_path,
      pattern = "cnmf_run\\.gene_spectra_tpm.*\\.dt_0_20\\.txt$",
      full.names = TRUE
    )
    if (length(file) == 0)
      stop(paste("No CNMF gene_spectra_tpm file found in:", dir_path))
    
    return(file[1])
  }
  
  # --- Internal helper: get top-N genes per program ---
  get_top_genes <- function(cnmf_matrix, top_n) {
    top_matrix <- matrix(nrow = nrow(cnmf_matrix), ncol = top_n)
    rownames(top_matrix) <- rownames(cnmf_matrix)
    
    for (i in 1:nrow(cnmf_matrix)) {
      row_data <- as.numeric(cnmf_matrix[i, ])
      top_idx <- order(row_data, decreasing = TRUE)[1:top_n]
      top_genes <- colnames(cnmf_matrix)[top_idx]
      top_matrix[i, ] <- top_genes
    }
    
    df <- as.data.frame(t(top_matrix))
    colnames(df) <- paste0("Usage_", colnames(df))
    return(df)
  }
  file1 <- find_cnmf_file(dir1)
  file2 <- find_cnmf_file(dir2)

  cnmf1 <- read.table(file1, header = TRUE, fill = TRUE, row.names = 1)
  cnmf2 <- read.table(file2, header = TRUE, fill = TRUE, row.names = 1)

  top_df1 <- get_top_genes(cnmf1, top_n)
  top_df2 <- get_top_genes(cnmf2, top_n)

  if (vf_only) {
    vfs <- VariableFeatures(object)
    
    top_df1 <- top_df1 %>%
      dplyr::mutate(across(everything(), ~ ifelse(.x %in% vfs, .x, NA)))
    
    top_df2 <- top_df2 %>%
      dplyr::mutate(across(everything(), ~ ifelse(.x %in% vfs, .x, NA)))
  }

  k1 <- ncol(top_df1)
  k2 <- ncol(top_df2)
  
  overlap_matrix <- matrix(0, nrow = k1, ncol = k2)
  
  for (i in 1:k1) {
    genes1 <- as.character(unlist(top_df1[, i]))
    for (j in 1:k2) {
      genes2 <- as.character(unlist(top_df2[, j]))
      overlap_matrix[i, j] <- length(intersect(genes1, genes2))
    }
  }
  
  rownames(overlap_matrix) <- paste0(label1, "_", 1:k1)
  colnames(overlap_matrix) <- paste0(label2, "_", 1:k2)

  p <- pheatmap::pheatmap(
    overlap_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%s",
    angle_col = 45,
    main = heatmap_title
  )
  
  print(p)

  return(list(
    overlap_matrix = overlap_matrix,
    top_df1 = top_df1,
    top_df2 = top_df2
  ))
}

#' Plot Paired Gene Expression with Links and Faceting
#'
#' Generates boxplots of gene expression (pre vs. post timepoint) with
#' individual data points linked for paired samples. Supports faceting by cohort.
#'
#' @param object A Seurat object containing expression data and metadata.
#' @param gene Character string specifying the gene to plot.
#' @param timepoint_col Character string specifying the metadata column for timepoints (e.g., "Timepoint").
#' @param participant_col Character string specifying the metadata column for patient ID (e.g., "Patient").
#' @param pre_label Character string for the 'pre' timepoint label (e.g., "Baseline").
#' @param post_label Character string for the 'post' timepoint label (e.g., "Post-Treatment").
#' @param facet_col Character string specifying the metadata column to facet by (e.g., "Cohort"). Default is NULL.
#' @param value_type Character string indicating whether to plot the mean ("mean") or median ("median") 
#'                   expression per patient. Default is "mean".
#'
#' @return A ggplot object showing paired expression changes.
#' @export
boxplot_paired_points <- function(
    object,
    gene,
    timepoint_col = "Timepoint",
    participant_col = "Patient",
    pre_label = "Baseline",
    post_label = "Post-Treatment",
    facet_col = NULL,
    value_type = "mean"
) {
  # 1. Input Validation
  if (!gene %in% rownames(object)) {
    stop(paste0("Gene '", gene, "' not found in Seurat object."))
  }
  
  # 2. Extract Data
  expr_matrix <- Seurat::GetAssayData(object, layer = "data")
  gene_expr <- expr_matrix[gene, ]
  
  metadata <- object@meta.data
  metadata$Expression <- gene_expr
  metadata$Patient <- metadata[[participant_col]]
  metadata$Timepoint <- metadata[[timepoint_col]]
  
  # 3. Filter and Aggregate Data per Patient
  
  # Filter for pre/post timepoints only
  plot_data <- metadata %>%
    dplyr::filter(.data$Timepoint %in% c(pre_label, post_label)) %>%
    dplyr::group_by(.data$Patient, .data$Timepoint, 
                    !!!if(!is.null(facet_col)) list(rlang::sym(facet_col)) else list()) %>%
    dplyr::summarize(
      Value = if (value_type == "mean") mean(Expression) else median(Expression),
      .groups = 'drop'
    ) %>%
    dplyr::ungroup()
  
  # Ensure Timepoint is a factor with the correct order
  plot_data$Timepoint <- factor(plot_data$Timepoint, levels = c(pre_label, post_label))
  
  # 4. Identify Paired Patients
  paired_data <- plot_data %>%
    dplyr::group_by(.data$Patient) %>%
    dplyr::filter(n() == 2) %>%
    dplyr::ungroup()
  
  # 5. Create the Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Timepoint, y = .data$Value)) +
    
    # Add boxplots for overall distribution
    ggplot2::geom_boxplot(
      width = 0.5, 
      outlier.shape = NA, # Remove outliers drawn by boxplot
      alpha = 0.6
    ) +
    
    # Add paired lines (geom_line uses paired_data)
    ggplot2::geom_line(
      data = paired_data, 
      ggplot2::aes(group = .data$Patient), 
      color = "gray", 
      linewidth = 0.5
    ) +
    
    # Add individual data points (geom_point uses all plot_data)
    ggplot2::geom_point(
      size = 2,
      # Use color aesthetic from Timepoint for visual distinction
      ggplot2::aes(color = .data$Timepoint) 
    ) +
    
    # Customize aesthetics
    ggplot2::labs(
      title = gene,
      x = rlang::as_label(rlang::sym(timepoint_col)),
      y = paste(stringr::str_to_title(value_type), "Expression")
    ) +
    ggplot2::theme_minimal()
  
  # 6. Add Faceting if requested
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)))
  }
  
  return(list(plot = p,
              data = paired_data))
}