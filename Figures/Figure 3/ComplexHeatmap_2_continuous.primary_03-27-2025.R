ComplexHeatMap_2 <- function(object,
                             data.type = 'CellToCell',
                             primary = 'Pseudotime',
                             secondary = 'CellType_Prox_Subset', # Cell type annotations
                             primary.cols = NULL,
                             secondary.cols = NULL,
                             features = NULL,
                             labels = NULL,
                             selected.row.anotations = NULL,
                             selected.label.size = 8,
                             use.scale.data = TRUE,
                             range.frac = 1,
                             row.dendrogram = FALSE) {
  # Packages
  require(tidyverse)
  require(RColorBrewer)
  require(ggplot2)
  require(ComplexHeatmap)
  require(circlize)
  
  # Add randomization column for ordering
  object$random <- sample(ncol(object))
  
  # Extract metadata
  meta.data <- object@meta.data
  meta.data$barcode <- rownames(meta.data)
  meta.data <- meta.data[order(
    meta.data[[primary]],
    meta.data[[secondary]],
    meta.data[['random']]
  ), ]
  
  # Get data for heatmap
  if (use.scale.data) {
    to.plot <- as.matrix(object@assays[[data.type]]$scale.data[features, meta.data$barcode])
  } else {
    to.plot <- as.matrix(object@assays[[data.type]]$data[features, meta.data$barcode])
  }
  
  # Set row and column names
  rownames(to.plot) <- features
  colnames(to.plot) <- NULL
  
  # Colors for primary (continuous variable)
  if (is.null(primary.cols)) {
    primary.cols <- circlize::colorRamp2(
      seq(min(meta.data[[primary]]), max(meta.data[[primary]]), length.out = 10),
      c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", 
        "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")
    )
  }
  
  # Colors for secondary (categorical variable)
  if (is.null(secondary.cols)) {
    cols.2 <- RColorBrewer::brewer.pal(n = 9, name = 'Set1')
    secondary.colors <- colorRampPalette(cols.2)(length(unique(meta.data[[secondary]])))
    names(secondary.colors) <- unique(meta.data[[secondary]])
  } else {
    secondary.colors <- secondary.cols[unique(meta.data[[secondary]])]
    names(secondary.colors) <- unique(meta.data[[secondary]])
  }
  
  # Define annotations
  stuff <- data.frame(
    primary = meta.data[[primary]],
    secondary = meta.data[[secondary]],
    check.names = FALSE
  )
  names(stuff) <- labels
  
  # Define colors for annotations
  colors <- list(
    primary = primary.cols,
    secondary = secondary.colors
  )
  names(colors) <- labels
  
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = stuff,
    col = colors,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9),  # ← Change this to desired label size
    gp = gpar(fontsize = 9)  # ← Controls text size *inside* the bar (optional)
  )
  
  # Color function for heatmap
  if (use.scale.data) {
    col_fun <- circlize::colorRamp2(
      c(min(to.plot) * range.frac, 0, max(to.plot) * range.frac),
      c("#440154", "#21918c", "#fde725")
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(0, max(to.plot) * range.frac),
      c("#440154", "#fde725")
    )
  }
  
  # Legend title
  legend.title <- if (use.scale.data) "Scaled Expression" else "Expression"
  
  # Row annotations if selected
  if (!is.null(selected.row.anotations)) {
    HAleft <- rowAnnotation(foo = anno_mark(
      at = which(rownames(to.plot) %in% selected.row.anotations),
      side = 'left',
      labels = rownames(to.plot)[rownames(to.plot) %in% selected.row.anotations],
      labels_gp = gpar(fontsize = selected.label.size)
    ))
  } else {
    HAleft <- NULL
  }
  
  # Create heatmap
  heatmap.object <- Heatmap(
    to.plot,
    column_split = stuff[, 1],
    column_gap = unit(0.3, 'mm'),
    col = col_fun,
    use_raster = FALSE,
    cluster_rows = row.dendrogram,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_row_dend = row.dendrogram,
    top_annotation = column_ha,
    left_annotation = HAleft,
    name = legend.title,
    row_names_side = "left",  # Labels only on the left
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(title_position = 'leftcenter-rot'),
    row_title = 'Differentially Expressed Genes (DEGs)',
    column_title = NULL
  )
  
  # Draw heatmap with legends right
  ht <- draw(
    heatmap.object,
    heatmap_legend_side = 'right',
    annotation_legend_side = 'right'
  )
  
  row.order.output <<- row_order(ht)
}
