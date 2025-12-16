#' Filter Genes by Minimum Cell Expression
#'
#' This function filters genes in a Seurat object based on the number of cells
#' in which they are expressed. Only genes expressed in at least `gene_filter`
#' cells are retained.
#'
#' @param object A Seurat object.
#' @param gene_filter Integer. Minimum number of cells a gene must be expressed in to be kept. Default is 50.
#' @param assay Character. Name of the assay to use for filtering. Default is "RNA".
#'
#' @return A Seurat object containing only the genes that meet the expression threshold.
#' @export
filter_genes_by_expression <- function(object, gene_filter = 50, assay = "RNA") {
  
  # Get raw counts
  # gene_expression_matrix <- object@assays[[assay]]$counts
  gene_expression_matrix <- Seurat::GetAssayData(object, assay = assay, layer = "counts")

  # Count cells with non-zero expression per gene
  gene_expression_count <- Matrix::rowSums(gene_expression_matrix > 0)

  # Keep genes expressed in at least `gene_filter` cells
  genes_to_keep <- base::names(gene_expression_count[gene_expression_count >= gene_filter])

  # Subset object to keep only selected genes
  object <- subset(object, features = genes_to_keep)
  
  message("Gene filtering done.")
  
  return(object)
}



#' Find integration anchors and set variable features
#'
#' This function aims to use anchor features as variable features for the object.
#'
#' @param object A Seurat object.
#' @param split.by Character string specifying the metadata column to split by.
#'   Default is "Patient"
#' @param selection.method Variable feature selection method for individual object. Default is "vst".
#' @param nfeatures Number of features to use as integration anchors. Default
#'   uses the function parameter nfeatures
#'
#' @return A Seurat object with anchor features set as variable features
#'
#' @examples
#' \dontrun{
#' seurat_obj <- anchor_features(seurat_obj, split.by = "Patient", nfeatures = 2000)
#' }
#'
#' @export
anchor_features <- function(object, split.by = "Patient", nfeatures = 2000, selection.method = "vst"){
  
  object <- NormalizeData(object)
  patient.list <- SplitObject(object, split.by = split.by)
  
  for (i in 1:length(patient.list)) {
    patient.list[[i]] <- NormalizeData(patient.list[[i]], verbose = FALSE)
    patient.list[[i]] <- FindVariableFeatures(patient.list[[i]], selection.method = selection.method, nfeatures = 2000, verbose = FALSE)
  }
  patient.list <- patient.list[sapply(patient.list, function(x) ncol(x) > 30)]
  anchors <- FindIntegrationAnchors(object.list = patient.list, anchor.features = nfeatures)
  VariableFeatures(object) <- anchors@anchor.features
  
  message("Define Variable Features by anchor features done.")
  
  return(object)
}

#' Preprocess Seurat Object.
#'
#' This function normalizes, selects variable features, and scales a Seurat object.
#' based on the chosen variable feature method and scaling method.
#'
#' @param object A Seurat object.
#' @param vf.method Character. Method for variable feature selection. Options: "vst", "dispersion", "mean.var.plot", "pseudobulk", "sct", "anchor".
#' @param nfeatures Integer. Number of variable features to select. Default is 2000.
#' @param regress.out Character vector. Variables to regress out during scaling. Default is NULL.
#' @param model.use Character. Model to use in scaling ("linear" or "poisson"). Default is "linear".
#' @param pseudobulk_var Character. Identity variable to use for pseudobulk aggregation. Required if vf.method = "pseudobulk".
#' @param pseudo_cutoff Numeric. Minimum cells per identity to include in pseudobulk. Default is 30.
#'
#' @return A Seurat object processed according to the chosen variable feature method.
#' @export
preprocess_obj <- function(object,
                                      vf.method = "vst",
                                      nfeatures = 2000,
                                      regress.out = NULL,
                                      model.use = "linear",
                                      pseudobulk_var = NULL,
                                      pseudo_cutoff = 30) {

  if(vf.method %in% c("vst", "dispersion", "mean.var.plot")) {
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, selection.method = vf.method, nfeatures = nfeatures)
    object <- ScaleData(object, vars.to.regress = regress.out, model.use = model.use)

  } else if(vf.method == "pseudobulk") {
    if(is.null(pseudobulk_var)) stop("'pseudobulk_var' must be provided for vf.method = 'pseudobulk'")
    
    Idents(object) <- pseudobulk_var
    ident_counts <- table(Idents(object))
    valid_idents <- names(ident_counts[ident_counts > pseudo_cutoff])
    
    if (length(valid_idents) == 0) {
      stop("No identities have more than ", pseudo_cutoff, " cells. ",
           "Try lowering pseudo_cutoff or use a different vf.method.")
    }
    obj2 <- subset(object, idents = valid_idents)

    avg_expr <- AverageExpression(obj2, assays = "RNA", layer = "counts", return.seurat = TRUE)
    avg_expr <- NormalizeData(avg_expr)
    avg_expr <- FindVariableFeatures(avg_expr, selection.method = "dispersion", nfeatures = nfeatures)

    object <- NormalizeData(object)
    VariableFeatures(object) <- VariableFeatures(avg_expr)
    object <- ScaleData(object, model.use = model.use)

  } else if(vf.method == "sct") {
    object <- SCTransform(object, variable.features.n = nfeatures, vars.to.regress = regress.out)

  } else {
    object <- anchor_features(object, nfeatures = nfeatures)
    object <- ScaleData(object, vars.to.regress = regress.out, model.use = model.use)
  }
  
  message("Preprocess object done.")

  return(object)
}


#' Select principal components based on variance explained
#'
#' This function analyzes the variance explained by each principal component
#' and selects the optimal number of PCs based on changes in variance explained.
#' It also generates an elbow plot with a red line indicating the selected PCs.
#'
#' @param object A Seurat object.
#' @param improved_diff_quantile Quantile threshold for changes in variance
#' explained. Lower values select more PCs. Default is 0.6.
#'
#' @return A numeric vector of selected PC indices
#' @examples
#' \dontrun{
#' selected_pcs <- select_PCs(seurat_obj, improved_diff_quantile = 0.6)
#' }
#'
#' @export
select_PCs <- function(object, improved_diff_quantile = 0.6){
  eigValues <- (object@reductions$pca@stdev)^2
  varExplained <- eigValues / sum(eigValues)
  dims <- 1:(floor(max(which(diff(varExplained) < quantile(diff(varExplained), 1-improved_diff_quantile)) + 1)))
  
  p <- ElbowPlot(object, ndims=ncol(object@reductions$pca@cell.embeddings)) +
    geom_vline(xintercept = length(dims), color="red")
  
  png("elbow.png")
  print(p)
  dev.off()
  message(paste0("Elbow plot saved to ", getwd(), "/elbow.png"))
  
  print(p)
  
  return(dims)
}

#' Plot mean-variance relationship with feature highlighting
#'
#' Generates mean-variance plots highlighting different sets of features:
#' variable features, selected features, and their overlap. Creates multiple
#' PNG files with the specified prefix.
#'
#' @param object A Seurat object.
#' @param features Character vector of selected feature names to
#'   highlight. If NULL, only the variable features plot is generated.
#'   Default is NULL
#' @param file_name Character string file_name for output PNG filenames.
#'   Default is "Mean-variance"
#'
#' @return Invisibly returns a data frame with gene statistics (mean, variance,
#'   log-transformed values)
#'
#' @examples
#' \dontrun{
#' # Plot only variable features
#' plot_mean_variance(seurat_obj)
#'
#' # Plot with selected features
#' selected_genes <- c("Gene1", "Gene2", "Gene3")
#' gene_stats <- plot_mean_variance(seurat_obj, selected_genes, file_name = "my_analysis")
#' }
#' @export
plot_mean_variance <- function(object, features = NULL, 
                               file_name = "Mean-variance") {
  # Compute mean and variance
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  gene_means <- Matrix::rowMeans(counts)
  gene_vars  <- Matrix::rowMeans(counts^2) - gene_means^2
  gene_vars[gene_vars < 0] <- 0

  gene_stats <- data.frame(
    gene = rownames(counts),
    mean = gene_means,
    variance = gene_vars,
    log_mean = log10(gene_means + 1),
    log_variance = log10(gene_vars + 1)
  )

  # Helper function to generate each plot
  make_plot <- function(gene_stats, highlight_genes, filename, title, label_genes = FALSE) {
    gene_stats$highlight <- ifelse(gene_stats$gene %in% highlight_genes, "Yes", "No")
    gene_stats <- gene_stats[order(gene_stats$highlight), ]

    p <- ggplot(gene_stats, aes(x = log_mean, y = log_variance)) +
      geom_point(aes(color = highlight), alpha = 0.4, size = 0.5) +
      scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
      theme_minimal(base_size = 12) +
      labs(
        title = title,
        x = "Log10 Mean Expression",
        y = "Log10 Variance",
        color = "Highlighted Genes"
      )

    if (label_genes) {
      p <- p + ggrepel::geom_text_repel(
        data = subset(gene_stats, highlight == "Yes"),
        aes(label = gene),
        size = 3, max.overlaps = 100,
        box.padding = 0.3, point.padding = 0.3,
        segment.color = "grey50"
      )
    }

    png(filename, width = 1600, height = 1600, res = 200)
    print(p)
    dev.off()
    
    return(p)
  }

  # 1. Variable features
  make_plot(
    gene_stats = gene_stats,
    highlight_genes = VariableFeatures(object),
    filename = sprintf("%s_variable.png", file_name),
    title = "Mean-Variance Plot with Variable Features"
  )

  # 2. Selected features
  if (!is.null(features)) {
    # Check which features exist
    existing_features <- intersect(features, rownames(counts))
    missing_features <- setdiff(features, rownames(counts))
    
    if (length(missing_features) > 0) {
      warning("The following selected features were not found in the object: ",
              paste(missing_features, collapse = ", "))
    }
    
    if (length(existing_features) == 0) {
      warning("None of the selected features were found in the object. Skipping selected feature plots.")
    } else {
      make_plot(gene_stats = gene_stats,
        highlight_genes = existing_features,
        filename = sprintf("%s_selected.png", file_name),
        title = "Mean-Variance Plot with Selected Features",
        label_genes = TRUE
      )
      # 3. Intersection
      overlap_genes <- intersect(VariableFeatures(object), existing_features)
      
      if (length(overlap_genes) > 0) {
        make_plot(
          gene_stats = gene_stats,
          highlight_genes = overlap_genes,
          filename = sprintf("%s_overlap.png", file_name),
          title = "Mean-Variance Plot with Overlapping Features",
          label_genes = TRUE
        )
      } else {
        message("No overlap between variable features and selected features.")
      }
    }
  }
  
  message("Plot Mean-Variance done.")

  invisible(gene_stats)
}


#' Export top genes per principal component
#'
#' @param object A Seurat object.
#' @param reduction Character; name of the dimensional reduction (default "pca")
#' @param top_n Integer; number of top positive/negative genes per PC (default 20)
#' @param n_pcs Integer; number of PCs to consider (default 30)
#' @param file_name Character; csv file to save results (default "PCs_top-genes.csv")
#' @return A data.frame of top genes per PC (invisible)
#' @export
export_top_pc_genes <- function(object,
                                reduction = "pca",
                                top_n = 20,
                                n_pcs = 30,
                                file_name = "PCs_top-genes.csv") {

  if (!reduction %in% names(object@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }

  # Extract loadings
  pca_loadings <- Loadings(object, reduction = reduction)

  # Check if requested n_pcs exceeds available PCs
  n_pcs <- min(n_pcs, ncol(pca_loadings))

  # Collect top genes for each PC
  top_genes_df <- do.call(rbind, lapply(1:n_pcs, function(pc) {
    loadings_pc <- pca_loadings[, pc]
    pos_genes <- names(sort(loadings_pc, decreasing = TRUE))[1:top_n]
    neg_genes <- names(sort(loadings_pc, decreasing = FALSE))[1:top_n]
    data.frame(
      PC = paste0("PC_", pc),
      Positive_Genes = pos_genes,
      Negative_Genes = neg_genes,
      stringsAsFactors = FALSE
    )
  }))

  # Save to CSV
  write.csv(top_genes_df, file_name, row.names = FALSE)
  message("Top genes per PC exported to ", file_name)
  head(top_genes_df)
}



#' Select marker genes based on product score
#'
#' This function selects marker genes from a differential expression results table
#' using a scoring system that combines statistical significance (`p_val_adj`) and
#' effect size (`avg_log2FC`). It can take results from both `FindAllMarkers` (> 2 groups) and `FindMarkers` (2 groups).
#' 
#' If using results from `FindMarkers`, `direction` doesn't take effect. Please make sure you add `cluster` column indicating the group name in the input.
#'
#' @param markers A data frame or tibble containing marker gene statistics.
#'   Must include the following columns:
#'   \describe{
#'     \item{gene}{Gene name or identifier.}
#'     \item{cluster}{Cluster or group label.}
#'     \item{p_val_adj}{Adjusted p-value (typically from a Seurat differential expression test).}
#'     \item{avg_log2FC}{Average log2 fold change between the cluster and others.}
#'   }
#' @param top_n Integer specifying the number of top marker genes to select per cluster,
#'   ranked by score. Default is 25.
#' @param direction Character string specifying whether to select upregulated
#'   (`"up"`) or downregulated (`"down"`) genes. Default is `"up"`.
#' @param adj_p_cutoff Numeric value specifying the statistical significance cutoff.
#' @param log2fc_cutoff Numeric value specifying the minimum absolute log2 fold change
#'   threshold. Default is 0.
#'
#' @details
#' The function computes a composite score for each gene as:
#' \deqn{score = -log10(p\_val\_adj) * |avg\_log2FC|}
#'
#' It then filters genes based on the specified direction (`avg_log2FC > log2fc_cutoff` for `"up"`,
#' `avg_log2FC < -log2fc_cutoff` for `"down"`), selects the top `n` genes by score per cluster. If `p_val_adj` is zero, it is replaced with
#' the smallest positive representable value to avoid infinite scores.
#'
#' @return A tibble containing the selected genes with columns:
#' \describe{
#'   \item{gene}{Gene name.}
#'   \item{cluster}{Cluster identifier.}
#'   \item{p_val_adj}{Adjusted p-value.}
#'   \item{avg_log2FC}{Average log2 fold change.}
#'   \item{score}{Computed score for ranking.}
#' }
#'
#' @examples
#' \dontrun{
#' markers <- FindAllMarkers(seurat_object)
#' selected_genes <- c("CD3D", "MS4A1", "LYZ")
#'
#' # Select top 20 upregulated genes with log2FC > 0.5
#' selected <- select_marker_genes_score(
#'   markers,
#'   top_n = 20,
#'   direction = "up",
#'   log2fc_cutoff = 0.5
#' )
#' head(selected)
#' }
#'
#' @export
select_marker_genes_score <- function(markers,
                                      top_n = 25,
                                      direction = "up",
                                      adj_p_cutoff = 0.05,
                                      log2fc_cutoff = 0) {
  markers$p_val_adj <- ifelse(markers$p_val_adj == 0, .Machine$double.xmin, markers$p_val_adj)
  markers$score <- (-log10(markers$p_val_adj)) * abs(markers$avg_log2FC)
  
  if(length(unique(markers$cluster)) > 2){
    if (direction == "up") {
      features <- markers %>%
        filter(p_val_adj <= adj_p_cutoff, avg_log2FC > log2fc_cutoff)
    } else {
      features <- markers %>%
        filter(p_val_adj <= adj_p_cutoff, avg_log2FC < -log2fc_cutoff)
    }
  }else{
    features <- markers %>%
      filter(p_val_adj <= adj_p_cutoff)
  }
  

  top_features <- features %>%
    group_by(cluster) %>%
    top_n(top_n, score) %>%
    ungroup() %>%
    distinct(gene, cluster, .keep_all = TRUE) %>%
    arrange(cluster, desc(score)) %>%
    as.data.frame()
  
  return(top_features)
}


#' Generate and save a Seurat DotPlot dynamically
#'
#' @param object A Seurat object.
#' @param features A gene list.
#' @param idents The identity column to group cells. E.g. seurat_clusters.
#' @param file_name Output file name (e.g., "dotplot.png").
#' @param width_factor Width multiplier per gene (default: 25).
#' @param height_factor Height multiplier per cluster (default: 40).
#' @param res Plot resolution in dpi (default: 100).
#' @param low_color,mid_color,high_color Colors for the expression gradient.
#'
#' @return Saves a PNG file of the dot plot.
#' @export
dotplots_png <- function(
    object,
    features,
    idents,
    file_name,
    width_factor = 25,
    height_factor = 42,
    res = 100,
    low_color = "#0024d6",
    mid_color = "#b4b6bf",
    high_color = "#d91111",
    ...
) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (!idents %in% colnames(object@meta.data)) {
    stop("`idents` must be a column in meta data.")
  }
  
  Idents(object) <- idents
  groups <- unique(object[[idents]] %>% unlist() %>% as.character())
  features <- unique(features)
  
  p <- DotPlot(object, features = features) +
    RotatedAxis() +
    scale_colour_gradient2(low = low_color, mid = mid_color, high = high_color)

  png(
    file_name,
    width = length(features) * width_factor,
    height = length(groups) * height_factor,
    res = res
  )
  print(p)
  dev.off()
  
  return(p)
}


#' Generate DotPlot PDFs for a Seurat Object
#'
#' This function creates a PDF containing two DotPlots for a given Seurat object:
#' one using a set of selected features, and another using top marker genes identified
#' by `FindAllMarkers()` and selected with `select_marker_genes_score()`.
#'
#' @param object A Seurat object.
#' @param features.1 A character vector of feature names to include in the first DotPlot.
#' @param features.2 A character vector of feature names to include in the second DotPlot.
#' @param file_name A string giving the output pdf filename (`dotplot.pdf`).
#'
#' @return A list containing:
#' \describe{
#'   \item{markers}{The full marker table from FindAllMarkers()}
#'   \item{genes}{The selected marker genes returned by select_marker_genes_score()}
#' }
#'
#' @examples
#' \dontrun{
#' dotplots_pdf(object = seurat_obj,
#'                      features.1 = c("GeneA", "GeneB"),
#'                      features.2 = c("GeneC", "GeneD"),
#'                      file_name = "dotplot.pdf")
#' }
#' @export
dotplots_pdf <- function(object, features.1, features.2, file_name = "dotplot.pdf") {
  
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object")
  
  if (!is.character(features.1) || length(features.1) == 0)
    stop("feature.1 must be a non-empty character vector")
  
  if (!is.character(features.2) || length(features.2) == 0)
    stop("feature.2 must be a non-empty character vector")
  
  pdf(
    file_name,
    width = length(unique(Idents(object))) * 5.7,
    height = length(unique(Idents(object))) * 0.6
  )

  # First dot plot
  print(
    DotPlot(object, features = features.1) +
      RotatedAxis() +
      scale_colour_gradient2(low = "#0024d6", mid = "#b4b6bf", high = "#d91111")
  )

  # Second dot plot
  print(
    DotPlot(object, features = unique(features.2)) +
      RotatedAxis() +
      scale_colour_gradient2(low = "#0024d6", mid = "#b4b6bf", high = "#d91111")
  )

  dev.off()

  # head(genes)
  # return(list(markers = markers, genes = genes))
}



#' Remove low-quality clusters based on feature count thresholds
#'
#' @param object A Seurat object.
#' @param lower_bound_param Numeric; the number of standard deviations below the mean
#'   nFeature_RNA that a cluster must fall to be considered low-quality. To be more stringent, decrease this parameter.
#'   Default is 1.
#' @param higher_bound_param Numeric; the number of standard deviations above the mean
#'   nFeature_RNA that a cluster must exceed to be considered low-quality. To be more stringent, decrease this parameter.
#'   Default is 2.
#'
#' @return A Seurat object with low-quality clusters removed.
#' @export
remove_low_quality_clusters <- function(object,
                                        lower_bound_param = 1,
                                        higher_bound_param = 2) {
  # --- Input validation ---
  if (!inherits(object, "Seurat")) stop("Input must be a Seurat object.")
  if (!all(c("seurat_clusters", "nFeature_RNA") %in% colnames(object@meta.data))) {
    stop("Metadata must contain 'seurat_clusters' and 'nFeature_RNA'.")
  }

  # --- Compute thresholds ---
  mean_nFeature <- mean(object$nFeature_RNA)
  sd_nFeature <- sd(object$nFeature_RNA)
  lower_bound <- mean_nFeature - lower_bound_param * sd_nFeature
  upper_bound <- mean_nFeature + higher_bound_param * sd_nFeature

  # --- Identify outlier clusters ---
  cluster_stats <- aggregate(nFeature_RNA ~ seurat_clusters, object@meta.data, mean)
  low_quality_clusters <- cluster_stats$seurat_clusters[
    cluster_stats$nFeature_RNA < lower_bound | cluster_stats$nFeature_RNA > upper_bound
  ]

  if (length(low_quality_clusters) == 0) {
    message("✅ No low-quality clusters detected. Returning all cells.")
  } else {
    message("Removing clusters: ", paste(low_quality_clusters, collapse = ", "))
  }

  # --- Visualization ---
  cluster_ids <- levels(Idents(object))
  cluster_colors <- setNames(rep("gray", length(cluster_ids)), cluster_ids)
  cluster_colors[as.character(low_quality_clusters)] <- "red"

  print(
    RidgePlot(object, features = "nFeature_RNA") +
      scale_fill_manual(values = cluster_colors) +
      ggtitle("Clusters flagged for removal (in red)")
  )

  # --- Subset Seurat object ---
  object <- subset(object, seurat_clusters %in% low_quality_clusters, invert = TRUE)
  
  return(object)
}


#' Automated parameter sweep for Seurat clustering
#'
#' This function performs iterative clustering across a range of neighbor
#' parameters (`k.param`) and resolution values. For each valid clustering result
#' (based on the number of clusters), it generates a `DimPlot` annotated with
#' cluster counts and saves all plots to a PDF file.
#'
#' @param object A Seurat object.
#' @param reduction dimensional reduction to be used. Default is "pca".
#' @param k_range Numeric vector specifying the range of `k.param` values to test
#'   in `FindNeighbors()`.
#' @param resolution_range Numeric vector specifying the range of resolution values
#'   to test in `FindClusters()`.
#' @param dims Integer vector of principal components to use in neighbor finding.
#' @param clusters_min Minimum acceptable number of clusters for a valid result.
#'   Defaults to a user-specified value.
#' @param clusters_max Maximum acceptable number of clusters for a valid result.
#'   Defaults to a user-specified value.
#' @param colors Optional vector of colors for clusters in plots.
#' @param file_name PDF file name to save all DimPlots. Default: "dimplot.pdf".
#'
#' @details
#' For each combination of `k.param` and `resolution`, the function runs:
#' \enumerate{
#'   \item \code{FindNeighbors(object, dims = dims, k.param = k)}
#'   \item \code{FindClusters(object, resolution = r)}
#' }
#' It then checks if the resulting number of clusters is between
#' `clusters_min` and `clusters_max`. If so, a `DimPlot` is created with
#' labeled clusters and saved to a combined PDF file named \file{dimplot.pdf}.
#'
#' Each accepted clustering is also stored in the Seurat object metadata
#' under a column named \code{"k<k>_r<r>"} (e.g., \code{"k20_r0.8"}).
#'
#' @return
#' A Seurat object with additional metadata columns for each valid clustering.
#' A PDF file (\file{dimplot.pdf}) containing all `DimPlot`s is also saved.
#' @examples
#' \dontrun{
#' seurat_obj <- clustering(
#'   object = seurat_obj,
#'   k_range = 10:30,
#'   resolution_range = seq(0.2, 1.2, by = 0.2),
#'   dims = 1:20,
#'   clusters_min = 5,
#'   clusters_max = 20
#' )
#' }
#'
#' @export
clustering <- function(
    object,
    reduction = "pca",
    k_range,
    resolution_range,
    dims,
    clusters_min = 5,
    clusters_max = 8,
    colors = NULL,
    file_name = "dimplot.pdf"
) {
  if (!inherits(object, "Seurat")) stop("Input must be a Seurat object.")
  if (missing(k_range) || missing(resolution_range) || missing(dims)) {
    stop("Please provide `k_range`, `resolution_range`, and `dims`.")
  }

  pdf(file_name, width = 4, height = 4)
  on.exit(dev.off(), add = TRUE)

  for (k in k_range) {
    message("---- Testing k.param = ", k, " ----")
    object <- FindNeighbors(object, dims = dims, k.param = k, reduction = reduction)

    for (r in resolution_range) {
      message("Running: k.param = ", k, " | resolution = ", r)
      object <- FindClusters(object, resolution = r)
      n_clusters <- length(unique(object$seurat_clusters))

      if (n_clusters >= clusters_min && n_clusters <= clusters_max) {
        title_text <- paste0("k.param: ", k, "  resolution: ", r)
        clust_counts <- table(object$seurat_clusters)
        clust_labels <- paste0("cluster", names(clust_counts), ": ", clust_counts)

        # Split long labels into multiple lines (5 per line)
        clust_str <- paste(split(clust_labels, ceiling(seq_along(clust_labels)/5)) |>
                             sapply(paste, collapse = "  "),
                           collapse = "\n")

        p <- DimPlot(object, label = TRUE, label.box = TRUE, cols = colors, alpha = 0.8) +
          labs(title = title_text, subtitle = clust_str) +
          theme(
            plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 7)
          )

        print(p)
        message("✅ Done: k.param = ", k, " | resolution = ", r)

        # Store cluster results
        temp <- paste0("k", k, "_r", r)
        object[[temp]] <- object$seurat_clusters
      } else if (n_clusters > clusters_max) {
        message("Too many clusters (", n_clusters, "); skipping higher resolutions.")
        next
      } else {
        message("Too few clusters (", n_clusters, "); skipping.")
        next
      }
    }
  }

  message("All clustering tests completed. Results saved to ", file_name)
  return(object)
}


#' Generate Clustered Heatmap for Seurat Object
#'
#' Scales data for selected genes, arranges cells by cluster, and produces a heatmap
#' with custom colors for clusters. Saves the heatmap as a PNG file.
#'
#' @param object A Seurat object.
#' @param features A character vector or a data frame with a column
#'   \code{gene} containing gene names to plot.
#' @param idents Character. Column name in \code{meta.data} containing identity.
#' @param colors Character vector of colors to use for clusters.
#' @param split_by_cluster Logical. If \code{TRUE}, splits the heatmap by cluster (default = TRUE).
#' @param cluster_rows Logical. Whether to cluster rows (genes) in the heatmap (default = FALSE).
#' @param cluster_columns Logical. Whether to cluster columns (cells) in the heatmap (default = FALSE).
#' @param file_name Character. File name for the saved heatmap PNG.
#'   If \code{NULL}, defaults to \code{"<cluster>_heatmap.png"}.
#' @param zlim Numeric vector of length 2 specifying limits for Z-score color scale (default = c(-2, 2)).
#' @param width Numeric. Width of the output PNG in inches (default = 6).
#' @param height Numeric. Height of the output PNG in inches (default = 8).
#' @param fontsize fontsize of rownames. Default is 3.
#' @param res Numeric. Resolution of the output PNG in DPI (default = 300).
#'
#' @return NULL. Heatmap is saved to file.
#' @export
heatmap_cell_level <- function(object,
                               features,
                               idents,
                               colors,
                               split_by_cluster = T,
                               cluster_rows = F,
                               cluster_columns = F,
                               file_name = NULL,
                               zlim = c(-2, 2),
                               width = 6,
                               height = NULL,
                               fontsize = 3,
                               res = 300) {

  meta <- object@meta.data
  meta[[idents]] <- factor(meta[[idents]], levels = sort(unique(meta[[idents]])))

  if (inherits(features, "character")) {
    genes_to_plot <- features
  }else{
    genes_to_plot <- features$gene
  }
  obj_hm <- ScaleData(object, features = genes_to_plot)

  mat <- GetAssayData(obj_hm, layer = "scale.data")
  mat <- mat[intersect(genes_to_plot, rownames(mat)), ]
  mat <- mat[, rownames(meta)]
  mat <- mat[rowSums(mat != 0) > 0, ]

  if(is.null(height)) {
    n_genes <- nrow(mat)
    height <- 8 + (n_genes - 120) * (2 / 50)
    height <- max(height, 6)
    message("Automatically calculated height: ", round(height, 2), " inches for ", n_genes, " genes")
  }

  group_colors <- colors[1:length(unique(meta[[idents]]))]
  names(group_colors) <- unique(meta[[idents]])

  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    Group = meta[[idents]],
    col = list(
      Group = group_colors
    ),
    annotation_name_side = "left"
  )

  col_fun <- circlize::colorRamp2(
    breaks = c(zlim[1], 0, zlim[2]),
    colors = c("blue", "white", "red")
  )

  if(split_by_cluster == T){
    
    p <- ComplexHeatmap::Heatmap(
      mat,
      name = "Expression",
      col = col_fun,
      top_annotation = col_anno,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_split = meta[[idents]],
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_title = "Genes",
      column_title = "Cells by groups",
      row_names_gp = grid::gpar(fontsize = fontsize),
      heatmap_legend_param = list(title = "Z-score")
    )
  }else{
    p <- ComplexHeatmap::Heatmap(
      mat,
      name = "Expression",
      col = col_fun,
      top_annotation = col_anno,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_split = NULL,
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_title = "Genes",
      column_title = "Cells by groups",
      row_names_gp = grid::gpar(fontsize = fontsize),
      heatmap_legend_param = list(title = "Z-score")
    )
  }
  
  if(!is.null(file_name)){
    png(file_name, width = width, height = height, units = "in", res = res)
    ComplexHeatmap::draw(p)
    dev.off()
    message("Heatmap saved to ", file_name)
  }
  
  return(p)
}

#' Plot Cluster Distribution Across Metadata Layers
#'
#' This function visualizes the distribution of Seurat clusters across one or more
#' metadata layers (e.g., sample, condition, or treatment). It supports one, two,
#' or three hierarchical layers. When three layers are provided, the function
#' generates a list of plots (one per level of the third layer).
#'
#' @param object A Seurat object. containing metadata columns specified in \code{idents}
#'   and \code{layers}.
#' @param idents A string specifying the metadata column that contains cluster assignments.
#' @param layers A character vector (length 1–3) specifying one or more metadata columns
#'   to group by when plotting. If three layers are provided, plots will be generated
#'   for each level of the third layer.
#' @param layer_orders An optional named list specifying the custom order of values
#'   for each metadata layer. The names must match entries in \code{layers}.
#' @param colors A vector of colors used to fill clusters in the plot.
#' @param percentage_by Optional character vector specifying the metadata column(s)
#'   by which percentages should be normalized. Defaults to \code{layers}.
#' @param title Optional plot title.
#' @param show_counts Logical indicating whether to display total cell counts on top of each bar.
#'   Default is FALSE.
#' @param count_size Numeric value for the size of count labels. Default is 3.
#' @param count_color Color for count labels. Default is "black".
#' @param count_vjust Vertical adjustment for count labels. Default is -0.5 (above bar).
#' @param .color_map Internal parameter for passing fixed color mapping in recursive calls.
#'
#' @return
#' \itemize{
#'   \item If \code{length(layers) <= 2}: A \code{ggplot} object showing cluster proportions.
#'   \item If \code{length(layers) == 3}: A named list of \code{ggplot} objects, one for each
#'     level of the third layer. The list has an attribute \code{"layer_order"} corresponding
#'     to the order of those levels.
#' }
#'
#' @details
#' The function calculates the percentage of cells per cluster within the specified
#' grouping layers. When three layers are provided, it subsets the Seurat object by
#' each level of the third layer and recursively generates plots using the first two
#' layers.
#'
#' @examples
#' \dontrun{
#' # Example: plot cluster distribution by sample and condition
#' stacked_bar_plots(
#'   object = seurat_obj,
#'   idents = "seurat_clusters",
#'   layers = c("sample", "condition"),
#'   layer_orders = list(
#'     sample = c("S1", "S2", "S3"),
#'     condition = c("Control", "Treatment")
#'   ),
#'   colors = RColorBrewer::brewer.pal(8, "Set2"),
#'   title = "Cluster distribution across samples",
#'   show_counts = TRUE
#' )
#' }
#'
#' @export
stacked_bar_plots <- function(object,
                                      idents,
                                      layers,
                                      layer_orders = list(),
                                      colors,
                                      percentage_by = NULL,
                                      title = NULL,
                                      show_counts = FALSE,
                                      count_size = 3,
                                      count_color = "black",
                                      count_vjust = -0.5,
                                      .color_map = NULL) {

  # --- Input validation ---
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!idents %in% colnames(object@meta.data)) {
    stop("Cluster column '", idents, "' not found in metadata")
  }

  if (length(layers) == 0 || length(layers) > 3) {
    stop("Number of layers must be between 1 and 3")
  }

  for (layer in layers) {
    if (!layer %in% colnames(object@meta.data)) {
      stop("Layer '", layer, "' not found in metadata")
    }
  }

  # --- Create fixed color map on first call ---
  if (is.null(.color_map)) {
    all_clusters <- unique(as.character(object[[idents, drop = TRUE]]))
    all_clusters_sorted <- stringr::str_sort(all_clusters, numeric = TRUE)
    .color_map <- colors[1:length(all_clusters_sorted)]
    names(.color_map) <- all_clusters_sorted
  }

  # --- Handle recursive case (3 layers) ---
  if (length(layers) == 3) {
    if(!is_empty(layer_orders) & layers[3] %in% names(layer_orders)){
      third_layer_values <- layer_orders[layers[3]][[1]]
    }else{
      third_layer_values <- unique(object[[layers[3]]] %>% unlist())
    }

    plot_list <- list()

    message("Input 3 layers. Use combine_stacked_bars() to combine them.")

    for (value in third_layer_values) {
      subobj <- subset(object, !!sym(layers[3]) == value)

      p <- stacked_bar_plots(
        object = subobj,
        idents = idents,
        layers = layers[1:2],
        layer_orders = layer_orders,
        colors = colors,
        percentage_by = percentage_by,
        title = paste(title, value, sep = " - "),
        show_counts = show_counts,
        count_size = count_size,
        count_color = count_color,
        count_vjust = count_vjust,
        .color_map = .color_map
      )

      plot_list[[as.character(value)]] <- p
    }

    attr(plot_list, "layer_order") <- third_layer_values
    return(plot_list)
  }

  # --- Main logic (1–2 layers) ---
  metadata <- object@meta.data

  # Calculate percentages
  group_vars <- c(layers, idents)
  summary_df <- metadata %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::tally() %>%
    dplyr::ungroup()

  if (is.null(percentage_by)) {
    percentage_grouping <- layers
  } else {
    percentage_grouping <- percentage_by
  }

  summary_df <- summary_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(percentage_grouping))) %>%
    dplyr::mutate(percentage = n / sum(n) * 100) %>%
    dplyr::ungroup()

  # Calculate total count per bar for labels
  if (show_counts) {
    total_counts <- summary_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(layers))) %>%
      dplyr::summarize(total_n = sum(n), .groups = "drop")

    summary_df <- summary_df %>%
      dplyr::left_join(total_counts, by = layers)
  }

  # Apply custom ordering to layers
  for (i in seq_along(layers)) {
    layer_name <- layers[i]

    if (layer_name %in% names(layer_orders)) {
      custom_levels <- layer_orders[[layer_name]]
      present_values <- unique(summary_df[[layer_name]])
      missing_values <- setdiff(present_values, custom_levels)

      if (length(missing_values) > 0) {
        warning(
          "Some values in layer '", layer_name, "' not in custom order: ",
          paste(missing_values, collapse = ", ")
        )
        custom_levels <- c(custom_levels, missing_values)
      }

      summary_df[[layer_name]] <- factor(summary_df[[layer_name]], levels = custom_levels)
    } else {
      summary_df[[layer_name]] <- factor(summary_df[[layer_name]])
    }
  }

  # Apply ordering to total_counts if needed
  if (show_counts) {
    for (i in seq_along(layers)) {
      layer_name <- layers[i]
      if (is.factor(summary_df[[layer_name]])) {
        total_counts[[layer_name]] <- factor(
          total_counts[[layer_name]],
          levels = levels(summary_df[[layer_name]])
        )
      }
    }
  }

  # Set cluster factor levels to ALL clusters (not just present ones)
  summary_df[[idents]] <- factor(
    summary_df[[idents]],
    levels = names(.color_map)
  )

  # Create plot
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = !!rlang::sym(layers[1]),
      y = percentage,
      fill = !!rlang::sym(idents)
    )
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(
      x = layers[1],
      y = "Percentage (%)",
      fill = "Group",
      title = title
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 70, hjust = 1),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black")
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    ggplot2::scale_fill_manual(
      values = scales::alpha(.color_map, 0.8),
      drop = FALSE
    )

  # Add total count labels on top of bars if requested
  if (show_counts) {
    p <- p + ggplot2::geom_text(
      data = total_counts,
      ggplot2::aes(
        x = !!rlang::sym(layers[1]),
        y = 100,
        label = total_n,
        fill = NULL
      ),
      size = count_size,
      color = count_color,
      vjust = count_vjust,
      inherit.aes = FALSE
    )
  }

  # Add facets if two layers
  if (length(layers) == 2) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", layers[2])), scales = "free_x", nrow = 1)
  }
  
  # print(p)
  return(p)
}

#' Combine a List of Cluster Distribution Plots
#'
#' This function takes a list of `ggplot` objects (typically generated by
#' `stacked_bar_plots()` for three-layer metadata) and combines them
#' into a single plot using `patchwork::wrap_plots()`. Legends are removed
#' from all plots except the last one, and subplot widths are dynamically
#' adjusted based on the number of bars in each plot.
#'
#' @param plot_list A named list of `ggplot` objects to be combined.
#' @param layer_order Optional character vector specifying the order of plots.
#'   If `NULL`, the function will attempt to use the `layer_order` attribute
#'   from `plot_list` if available.
#'
#' @return A combined `ggplot` object (patchwork object) with all subplots arranged
#'   horizontally and legends only on the last plot.
#'
#' @details
#' The function automatically:
#' \itemize{
#'   \item Removes legends from all but the last plot.
#'   \item Reorders plots according to `layer_order` or the `layer_order` attribute.
#'   \item Dynamically adjusts subplot widths based on the number of bars in each plot.
#' }
#' This is useful for visualizing cluster distributions across multiple subsets
#' of a dataset in a single figure.
#'
#' @examples
#' \dontrun{
#' # Combine a list of plots generated for each level of a third layer
#' plot_list <- stacked_bar_plots(
#'   object = seurat_obj,
#'   idents = "seurat_clusters",
#'   layers = c("sample", "condition", "treatment"),
#'   layer_orders = list(
#'     sample = c("S1", "S2"),
#'     condition = c("Control", "Treatment"),
#'     treatment = c("DrugA", "DrugB")
#'   ),
#'   colors = RColorBrewer::brewer.pal(8, "Set2")
#' )
#' combined <- combine_stacked_bars(plot_list)
#' print(combined)
#' }
#'
#' @export
combine_stacked_bars <- function(plot_list, layer_order = NULL) {
  if (length(plot_list) == 0) {
    stop("Plot list is empty")
  }

  # Use stored order from attributes if available and no explicit order provided
  if (is.null(layer_order) && !is.null(attr(plot_list, "layer_order"))) {
    layer_order <- attr(plot_list, "layer_order")
  }

  # Reorder plot_list if layer_order is provided
  if (!is.null(layer_order)) {
    valid_layers <- layer_order[layer_order %in% names(plot_list)]
    if (length(valid_layers) > 0) {
      plot_list <- plot_list[valid_layers]
    }
    missing_layers <- setdiff(names(plot_list), valid_layers)
    if (length(missing_layers) > 0) {
      plot_list <- c(plot_list[valid_layers], plot_list[missing_layers])
    }
  }

  # Remove legends from all plots except the last one
  if (length(plot_list) > 1) {
    plot_list[1:(length(plot_list) - 1)] <- lapply(plot_list[1:(length(plot_list) - 1)], function(p) {
      p + ggplot2::theme(legend.position = "none")
    })
  }

  # Calculate dynamic widths based on number of bars in each subplot
  widths <- sapply(plot_list, function(p) {
    plot_data <- ggplot2::ggplot_build(p)$data[[1]]
    if (!is.null(plot_data)) {
      n_bars <- sum(near(plot_data$y, 100, tol = 0.001), na.rm = TRUE)
      return(max(n_bars, 1))
    } else {
      return(1)
    }
  })

  # Normalize widths to sum to n_plots
  total_width <- sum(widths)
  if (total_width > 0) {
    widths <- (widths / total_width) * length(plot_list)
  }

  # Combine plots with dynamic widths
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = length(plot_list),
                                         widths = widths)

  return(combined_plot)
}

#' Optimize Single Cell RNA-seq Analysis
#'
#' This function performs a comprehensive single-cell RNA-seq analysis pipeline including
#' quality control, normalization, dimensionality reduction, clustering optimization, and
#' visualization. It automatically saves parameters, session info, and generates multiple
#' QC and analysis plots.
#'
#' @param object A Seurat object.
#' @param normalization.method Method for normalization. Default is "LogNormalize". Please refer to Seurat::NormalizeData() for more options.
#' @param vf.method Method for variable feature selection. Options include "vst", "dispersion",
#'   "pseudobulk", "sct", or "anchor". Default is "vst".
#' @param model.use Model to use for scaling. Default is "linear". Please refer to Seurat::ScaleData() for more options.
#' @param pseudobulk_var Variable name for pseudobulk aggregation. Default is "patient_tp_mutation".
#' @param pseudo_cutoff Minimum number of cells required for pseudobulk calculation. Default is 30.
#' @param colors Named vector of colors for plotting. Default is NULL.
#' @param mt.cutoff Maximum percent mitochondrial content threshold for cell filtering. Default is 25.
#' @param nfeatures Number of variable features to select. Default is 2000.
#' @param harmony Logical indicating whether to run Harmony batch correction. Default is FALSE.
#' @param harmony_vars Character vector of variable names to use for Harmony correction. e.g. run, experiment, patient. Default is NULL.
#' @param reduction Dimensionality reduction method to use ("pca" or "harmony"). Default is "pca".
#' @param improved_diff_quantile Quantile threshold for PC selection. Lower values select more PCs. Default is 0.6.
#' @param k_range Vector of k values for nearest neighbor calculations. Default is c(10, 20, 30).
#' @param resolution_range Vector of resolution values for clustering. Default is seq(0.1, 0.3, 0.1).
#' @param min.dist Minimum distance parameter for UMAP. Default is 0.4.
#' @param seed Random seed for reproducibility. Default is 17.
#' @param gene_filter Minimum number of cells a gene must be expressed in to be retained. Default is 50.
#' @param regress.out Variables to regress out during scaling, e.g. proliterating genes. Default is NULL.
#' @param clusters_min Minimum number of clusters expected. Default is 4.
#' @param clusters_max Maximum number of clusters expected. Default is 13.
#' @param min.pct Minimum percentage expressed required for a gene to become DEG. Default is 0.3.
#' @param lower_bound_param Lower bound parameter for low-quality cluster removal. To be more stringent, decrease this parameter. Default is 1.
#' @param higher_bound_param Upper bound parameter for low-quality cluster removal. To be more stringent, increase this parameter. Default is 2.
#' @param cluster_distribution_layers Character vector of metadata columns for cluster distribution plots. Default is NULL.
#' @param cluster_distribution_layer_orders List of factor level orders for distribution layers. Default is NULL.
#' @param out_dir Output directory for saving results. Default is current working directory.
#' @param features Data frame or list of selected features for visualization. Default is NULL.
#'
#' @return A processed Seurat object with clustering and dimensionality reduction completed.
#'   The object is also saved as "object.rds" in the output directory.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates output directory and logs parameters and session info
#'   \item Filters genes by expression threshold
#'   \item Processes variable features using specified method
#'   \item Runs PCA and optionally Harmony batch correction
#'   \item Performs initial clustering at high resolution for QC
#'   \item Removes low-quality clusters and high mitochondrial content cells
#'   \item Re-processes variable features after QC
#'   \item Selects optimal number of PCs
#'   \item Runs clustering across multiple k and resolution parameters
#'   \item Generates comprehensive visualizations including cluster distributions, dotplots, and heatmaps
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' seurat_obj <- optimize_single_cell(
#'   object = my_seurat_object,
#'   out_dir = "results/analysis_v1"
#' )
#'
#' # With Harmony batch correction
#' seurat_obj <- optimize_single_cell(
#'   object = my_seurat_object,
#'   harmony = TRUE,
#'   harmony_vars = c("batch", "donor"),
#'   out_dir = "results/harmony_analysis"
#' )
#'
#' # Custom clustering parameters
#' seurat_obj <- optimize_single_cell(
#'   object = my_seurat_object,
#'   k_range = c(15, 25, 35),
#'   resolution_range = seq(0.2, 0.8, 0.2),
#'   clusters_min = 5,
#'   clusters_max = 15,
#'   out_dir = "results/custom_clustering"
#' )
#' }
optimize_single_cell <- function(object = NULL,
                                 normalization.method = "LogNormalize",
                                 vf.method = "vst",
                                 model.use = "linear",
                                 pseudobulk_var = "patient_tp_mutation",
                                 pseudo_cutoff = 30,
                                 colors = NULL,
                                 mt.cutoff = 25,
                                 nfeatures = 2000,
                                 harmony = FALSE,
                                 harmony_vars = NULL,
                                 reduction = "pca",
                                 improved_diff_quantile = 0.6,
                                 k_range = c(10, 20, 30),
                                 resolution_range = seq(0.1, 0.3, 0.1),
                                 min.dist = 0.4,
                                 seed = 17,
                                 gene_filter = 50,
                                 regress.out = NULL,
                                 clusters_min = 4,
                                 clusters_max = 13,
                                 min.pct = 0.3,
                                 lower_bound_param = 1,
                                 higher_bound_param = 2,
                                 cluster_distribution_layers = NULL,
                                 cluster_distribution_layer_orders = NULL,
                                 out_dir = getwd(),
                                 features = NULL) {

  dir.create(out_dir)
  setwd(out_dir)

  # Capture all parameters (explicit + defaults)
  params <- as.list(match.call())[-1]  # Get call without function name
  params$object <- NULL                 # Remove large Seurat object

  # Get function defaults
  func_formals <- formals(sys.function(sys.parent()))
  func_formals$object <- NULL           # Remove object from defaults too
  # Merge: defaults first, then override with explicit calls
  all_params <- c(func_formals, params)
  all_params <- all_params[!duplicated(names(all_params), fromLast = TRUE)]

  # Evaluate to get actual values (not symbols/promises)
  all_params <- lapply(all_params, function(x) {
    if (is.symbol(x) || is.language(x)) {
      tryCatch(eval(x, envir = parent.frame(2)), error = function(e) x)
    } else {
      x
    }
  })

  # Remove any empty/NULL entries
  all_params <- all_params[!sapply(all_params, is.null)]
  
  safe_params <- all_params[sapply(all_params, function(x) {
    is.atomic(x) || is.list(x) || is.character(x) || is.numeric(x) || is.logical(x)
  })]

  # Save as YAML
  dir.create("logs", showWarnings = FALSE)
  param_file <- file.path(
    "logs",
    paste0("params_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".yaml")
  )
  yaml::write_yaml(safe_params, param_file)
  
  # --- Continue with your function code ---
  message("parameters saved to ", param_file)

  session_file <- file.path("logs", paste0("session_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  sink(session_file)
  print(sessioninfo::session_info())
  sink()
  message("Session info saved to ", session_file)

  pdf("qc.pdf")

  n1 <- dim(object)[1]
  object <- filter_genes_by_expression(object, gene_filter = gene_filter, assay = "RNA")
  n2 <- dim(object)[1]
  object <- preprocess_obj(object,
                                      vf.method = vf.method,
                                      nfeatures = nfeatures,
                                      regress.out = regress.out,
                                      model.use = model.use,
                                      pseudobulk_var = pseudobulk_var,
                                      pseudo_cutoff = pseudo_cutoff)
  object <- RunPCA(object, npcs = 30, verbose = F)
  if(harmony == T){
    object <- RunHarmony(object, group.by.vars = harmony_vars)
    reduction <- "harmony"
  }
  object <- RunUMAP(object, reduction = reduction, dims = 1:30, seed.use = seed)
  object <- FindNeighbors(object, dims = 1:30)
  object <- FindClusters(object, resolution = 2)

  print(DimPlot(object, label = TRUE, label.box = TRUE))
  print(FeaturePlot(object, features = "percent.mt"))
  print(FeaturePlot(object, features = "nCount_RNA"))

  n4 <- dim(object)[2]
  object <- remove_low_quality_clusters(object,
                                        lower_bound_param = lower_bound_param,
                                        higher_bound_param = higher_bound_param)

  object <- subset(object, percent.mt < mt.cutoff)

  n5 <- dim(object)[2]
  grid.newpage()
  grid.draw(grid.text(
    paste0("Cells before QC:", n4, "\n",
           "Cells after QC:", n5),
    x = 0.5, y = 0.6
  ))

  object <- filter_genes_by_expression(object, gene_filter = gene_filter)
  n3 <- dim(object)[1]
  grid.draw(grid.text(
    paste0("Genes before filtering:", n1, "\n",
           "Genes after 1st filtering:", n2, "\n",
           "Genes after 2nd filtering:", n3),
    x = 0.5, y = 0.4
  ))
  dev.off()

  object <- preprocess_obj(object,
                                      vf.method = vf.method,
                                      nfeatures = nfeatures,
                                      regress.out = regress.out,
                                      model.use = model.use,
                                      pseudobulk_var = pseudobulk_var,
                                      pseudo_cutoff = pseudo_cutoff)

  # export gene rank - conditionally applicable
  if(vf.method == "vst"){
    hvf <- HVFInfo(object)
    hvf$gene <- rownames(hvf)
    hvf <- hvf[order(hvf$variance.standardized, decreasing = T), ]
  }else if(vf.method == "dispersion") {
    hvf <- HVFInfo(object)
    hvf$gene <- rownames(hvf)
    hvf <- hvf[order(hvf$mvp.dispersion, decreasing = T), ]
  }else if(vf.method == "pseudobulk"){
    hvf <- HVFInfo(avg_expr)
    hvf$gene <- rownames(hvf)
    hvf <- hvf[order(hvf$mvp.dispersion, decreasing = T), ]
  }else if(vf.method == "sct"){
    hvf <- HVFInfo(object[["SCT"]])
  }
  hvf$rank <- 1:nrow(hvf)
  hvf$features <- "No"
  hvf$features[hvf$gene %in% features] <- "Yes"
  expr_matrix <- GetAssayData(object, layer = "counts")
  gene_cell_counts <- Matrix::rowSums(expr_matrix > 0)
  hvf$num_cells_expressed <- gene_cell_counts[match(hvf$gene, names(gene_cell_counts))]
  write.csv(hvf, "HVG_info.csv")

  plot_mean_variance(
    object = object,
    features = features
  )
  object <- RunPCA(object, npcs = 30)
  if(harmony == FALSE){
    export_top_pc_genes(object, reduction = reduction, top_n = 20, n_pcs = 30, file_name = "PCs_top-genes.csv")
  }

  dims <- select_PCs(object,
                     improved_diff_quantile = improved_diff_quantile)
  print(dims)

  if(harmony == T){
    object <- RunHarmony(object, group.by.vars = harmony_vars)
    reduction <- "harmony"
  }
  object <- RunUMAP(object, reduction = reduction, dims = dims, seed.use = seed)
  object <- clustering(object, reduction = reduction, colors = colors, k_range = k_range, resolution_range = resolution_range, dims = dims, clusters_min = clusters_min, clusters_max = clusters_max)
  clusters <- grep("^k\\d+_r\\d+\\.\\d+$", colnames(object@meta.data), value = TRUE)
  print(clusters)
  for(temp_clusters in clusters){
    Idents(object) <- temp_clusters
    print(table(Idents(object)))
    
    if(!is.null(cluster_distribution_layers)){
      plots <- stacked_bar_plots(
        object = object,
        idents = temp_clusters,
        layers = cluster_distribution_layers,
        layer_orders = cluster_distribution_layer_orders,
        colors = colors,
        title = "Cluster Distribution"
      )
      
      n_layers <- object[[cluster_distribution_layers]] %>% unlist() %>% as.character %>% unique() %>% length()
      width <- 5 + 0.3 * n_layers
      
      if(length(cluster_distribution_layers) == 3){
        combined_3layer_plot <- combine_stacked_bars(plots)
        ggsave(paste0(temp_clusters, "_barplot.png"), plot = combined_3layer_plot, width = width, height = 5.5, dpi = 300)
      }else{
        ggsave(paste0(temp_clusters, "_barplot.png"), plot = plots, width = width, height = 5.5, dpi = 300)
      }
    }
    
    markers <- FindAllMarkers(object, min.pct = min.pct)
    dge <- select_marker_genes_score(markers)

    dotplots_pdf(
      object = object,
      features.1 = features,
      features.2 = dge$gene,
      file_name = paste0(temp_clusters, "_dotplot.pdf")
    )

    heatmap_cell_level(object, 
                       features = dge$gene, 
                       idents = temp_clusters, 
                       colors, 
                       file_name = paste0(temp_clusters, "_heatmap.png"))
  }
  saveRDS(object, "object.rds")
  return(object)
}
