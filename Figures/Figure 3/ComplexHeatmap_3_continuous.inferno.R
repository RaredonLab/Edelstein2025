ComplexHeatMap_3.inferno <- function(object,
                             data.type = 'CellToCell',
                             primary = 'Pseudotime',
                             secondary = 'CellType.sub',
                             tertiary = NULL,
                             primary.cols = NULL,
                             secondary.cols = NULL,
                             tertiary.cols = NULL,
                             features = NULL,
                             labels = NULL,
                             selected.row.anotations = NULL,
                             selected.label.size = 10,
                             use.scale.data = TRUE,
                             range.frac = 1,
                             row.dendrogram = FALSE,
                             column_split = NULL) {
  
  require(tidyverse)
  require(RColorBrewer)
  require(ggplot2)
  require(ComplexHeatmap)
  require(circlize)
  require(viridis)  # For inferno
  
  object$random <- sample(ncol(object))
  meta.data <- object@meta.data
  meta.data$barcode <- rownames(meta.data)
  
  # Sort metadata by primary, secondary, tertiary (if any), and a random shuffle to stabilize ties
  meta.data <- meta.data[order(
    if (!is.null(primary)) meta.data[[primary]] else NULL,
    meta.data[[secondary]],
    if (!is.null(tertiary)) meta.data[[tertiary]] else NULL,
    meta.data[['random']]
  ), ]
  
  if (use.scale.data) {
    to.plot <- as.matrix(object@assays[[data.type]]$scale.data[features, meta.data$barcode])
  } else {
    to.plot <- as.matrix(object@assays[[data.type]]$data[features, meta.data$barcode])
  }
  
  rownames(to.plot) <- features
  colnames(to.plot) <- NULL
  
  # Primary (continuous) color scale
  if (is.null(primary.cols)) {
    primary.cols <- circlize::colorRamp2(
      seq(min(meta.data[[primary]]), max(meta.data[[primary]]), length.out = 10),
      rev(RColorBrewer::brewer.pal(10, "Spectral"))
    )
  }
  
  # Secondary (categorical)
  if (is.null(secondary.cols)) {
    cols.2 <- RColorBrewer::brewer.pal(9, 'Set1')
    secondary.colors <- colorRampPalette(cols.2)(length(unique(meta.data[[secondary]])))
    names(secondary.colors) <- unique(meta.data[[secondary]])
  } else {
    secondary.colors <- secondary.cols[unique(meta.data[[secondary]])]
    names(secondary.colors) <- unique(meta.data[[secondary]])
  }
  
  # Tertiary (categorical)
  if (!is.null(tertiary)) {
    if (is.null(tertiary.cols)) {
      cols.3 <- RColorBrewer::brewer.pal(8, 'Dark2')
      tertiary.colors <- colorRampPalette(cols.3)(length(unique(meta.data[[tertiary]])))
      names(tertiary.colors) <- unique(meta.data[[tertiary]])
    } else {
      tertiary.colors <- tertiary.cols[unique(meta.data[[tertiary]])]
      names(tertiary.colors) <- unique(meta.data[[tertiary]])
    }
  }
  
  # Build annotations
  stuff <- data.frame(
    primary = meta.data[[primary]],
    secondary = meta.data[[secondary]],
    check.names = FALSE
  )
  names(stuff) <- labels[1:2]
  
  colors <- list(
    primary = primary.cols,
    secondary = secondary.colors
  )
  names(colors) <- labels[1:2]
  
  if (!is.null(tertiary)) {
    stuff[[labels[3]]] <- meta.data[[tertiary]]
    colors[[labels[3]]] <- tertiary.colors
  }
  
  column_ha <- HeatmapAnnotation(
    df = stuff,
    col = colors,
    annotation_name_side = "left"
  )
  
  # ✅ Expression color scale using inferno
  inferno.colors <- viridis::inferno(256)
  col_fun <- if (use.scale.data) {
    circlize::colorRamp2(
      seq(min(to.plot) * range.frac, max(to.plot) * range.frac, length.out = 256),
      inferno.colors
    )
  } else {
    circlize::colorRamp2(
      seq(0, max(to.plot) * range.frac, length.out = 256),
      inferno.colors
    )
  }
  
  # Row annotation if selected
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
  
  # Set column split if not provided
  if (is.null(column_split)) {
    column_split <- NULL
  }
  
  heatmap.object <- Heatmap(
    to.plot,
    column_split = column_split,
    column_gap = unit(0.3, 'mm'),
    col = col_fun,
    use_raster = FALSE,
    cluster_rows = row.dendrogram,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    show_row_dend = row.dendrogram,
    top_annotation = column_ha,
    left_annotation = HAleft,
    name = if (use.scale.data) "Scaled Expression" else "Expression",
    row_names_side = "left",
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      title = if (use.scale.data) "Scaled Expression" else "Expression",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      title_position = 'leftcenter-rot'
    ),
    row_title = '',
    column_title = NULL
  )
  
  ht <- draw(
    heatmap.object,
    heatmap_legend_side = 'right',
    annotation_legend_side = 'right'
  )
  row.order.output <<- row_order(ht)
}
