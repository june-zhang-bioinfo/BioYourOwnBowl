# =============================================================================
# Package-wide imports (roxygen2 controlled)
# =============================================================================

# Seurat Functions
#' @importFrom Seurat DimPlot NormalizeData FindVariableFeatures FindIntegrationAnchors SCTransform ScaleData ElbowPlot FeaturePlot RidgePlot
#' @importFrom Seurat RunPCA RunUMAP FindNeighbors FindClusters SplitObject `Idents<-` `VariableFeatures<-`
#' @importFrom Seurat DefaultAssay VariableFeatures FindMarkers FindAllMarkers
#' @importFrom Seurat Idents GetAssayData AverageExpression Loadings RotatedAxis DotPlot HVFInfo FetchData

#' @importFrom scCustomize FeaturePlot_scCustom

# ggplot2
#' @importFrom ggplot2 ggsave geom_point geom_vline geom_hline geom_line geom_violin geom_dotplot geom_boxplot facet_wrap ggplot aes ggtitle
#' @importFrom ggplot2 scale_colour_gradient2 scale_color_manual scale_fill_manual labs theme theme_classic theme_minimal geom_text element_text element_rect

# dplyr
#' @importFrom dplyr %>% filter group_by summarize ungroup mutate top_n near desc across everything left_join select distinct arrange all_of

# tidyr/tibble (only if truly needed)
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble tibble rownames_to_column

# Other packages
#' @importFrom stringr str_sort
#' @importFrom circlize colorRamp2
#' @importFrom grid unit grid.newpage grid.draw grid.text gpar
#' @importFrom eulerr euler
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom grDevices colorRampPalette dev.list dev.off pdf png rainbow
#' @importFrom stats aggregate as.formula quantile sd setNames median
#' @importFrom utils write.csv read.table head combn
#' @importFrom glue glue
#' @importFrom ggpubr stat_compare_means
#' @importFrom jsonlite toJSON
#' @importFrom pheatmap pheatmap
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom anndata write_h5ad
#' @importFrom convert2anndata convert_to_anndata convert_seurat_to_sce
#' @importFrom ggrepel geom_text_repel
#' @importFrom matrixStats rowMedians rowMaxs rowMins
#' @importFrom sessioninfo session_info
#' @importFrom yaml write_yaml
#' @importFrom harmony RunHarmony
#' @importFrom rlang .data sym
#' @importFrom purrr is_empty

NULL

# =============================================================================
# Global variables
# =============================================================================

utils::globalVariables(c(
  "x", "y", "gene", "genes", "features", "avg_log2FC", "avg_log2FC_1", "avg_log2FC_2",
  "p_val_adj", "p_val_adj_1", "p_val_adj_2", "color_group", "score",
  "cluster", "cell_type", "mean_expression", "log_mean", "log_variance",
  "highlight", "percentage", "total_n", "patient_id", "sample_id",
  "cohort_tp", "Specificity", "participant", "seurat_clusters",
  "Column", "Value", "MinExpr", "MedianExpr", "MaxExpr", "percent.mt",
  "avg_expr", "py_code", "cell_id", "usage", "Expression", "n", "out_dir", "ident"
))
