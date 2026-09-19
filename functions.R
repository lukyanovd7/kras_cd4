#-----------------------------------------------------------------FUNCTIONS------------------------------------------------------


quality_plot_custom <- function(object, group.by, raster = FALSE, colors = NULL){
  
  mt_plot <- VlnPlot(object, features = c('prc.mt'), pt.size = 0, group.by=group.by, raster = raster, cols = colors)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  nCount_plot <- VlnPlot(object, features = c('nCount_RNA'), pt.size = 0, group.by=group.by, raster = raster, cols = colors)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  nFeature_plot <- VlnPlot(object, features = c('nFeature_RNA'), pt.size = 0, group.by=group.by, raster = raster, cols = colors)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    theme(legend.position='none')
  
  mt_plot/nCount_plot/nFeature_plot
}

plot_batch_fr <- function(seurat_obj, cluster_column, batch_column, palette, normalize = TRUE) {
  # Convert string column names to symbols
  cluster_sym <- rlang::sym(cluster_column)
  batch_sym <- rlang::sym(batch_column)
  
  # Count cells per cluster and batch
  data <- seurat_obj@meta.data %>%
    dplyr::select(!!cluster_sym, !!batch_sym) %>%
    dplyr::group_by(!!cluster_sym, !!batch_sym) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!batch_sym, values_from = count, values_fill = 0)
  
  # Normalize batch columns across clusters (optional)
  if (normalize) {
    batch_names <- setdiff(names(data), cluster_column)
    
    data_transformed <- data %>%
      dplyr::mutate(across(all_of(batch_names), ~ . / mean(.))) %>%
      dplyr::mutate(!!cluster_sym := data[[cluster_column]])
  } else {
    data_transformed <- data
  }
  
  # Calculate per-cluster fractions
  batch_frac <- data_transformed %>%
    rowwise() %>%
    mutate(rowSum = sum(c_across(-!!cluster_sym))) %>%
    ungroup() %>%
    mutate(across(-c(!!cluster_sym, rowSum), ~ . / rowSum)) %>%
    dplyr::select(-rowSum)
  
  # Reshape to long format
  plot_data <- batch_frac %>%
    tidyr::pivot_longer(cols = -!!cluster_sym, names_to = batch_column, values_to = "fraction")
  
  # Plot stacked bar chart
  ggplot(plot_data, aes(x = factor(!!cluster_sym), y = fraction, fill = !!rlang::sym(batch_column))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = cluster_column, y = "Fraction", fill = batch_column) +
    scale_y_continuous() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = palette)
}