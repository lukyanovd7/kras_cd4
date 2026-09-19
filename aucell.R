#-----------------------------------------------------------------LIBRARIES------------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(mltest)

#-----------------------------------------------------------------FUNCTIONS------------------------------------------------------


quality_plot_custom <- function(object, group.by){
  
  mt_plot <- VlnPlot(object, features = c('prc.mt'), pt.size = 0, group.by=group.by, raster = FALSE)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    xlab('')+
    theme(legend.position='none',
          axis.text.x = element_text(angle = 60))
  
  nCount_plot <- VlnPlot(object, features = c('nCount_RNA'), pt.size = 0, group.by=group.by, raster = FALSE)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    xlab('')+
    theme(legend.position='none',
          axis.text.x = element_blank())
  
  nFeature_plot <- VlnPlot(object, features = c('nFeature_RNA'), pt.size = 0, group.by=group.by, raster = FALSE)+
    stat_summary(fun = mean, geom='point', size = 3, colour = "red")+
    stat_summary(fun = median, geom='point', size = 3, colour = "orange")+
    xlab('')+
    theme(legend.position='none',
          axis.text.x = element_blank())
  
  #mt_plot|nCount_plot|nFeature_plot
  # ggpubr::ggarrange(nCount_plot + ggpubr::rremove('x.text'),
  #                   nFeature_plot + ggpubr::rremove('x.text'),
  #                   mt_plot,
  #                   ncol = 1)
  
  final_plot <- (nCount_plot / nFeature_plot / mt_plot)
  
  final_plot
}

plot_batch_fr <- function(seurat_obj, cluster_column, batch_column, palette, normalize = TRUE) {
  
  cluster_sym <- rlang::sym(cluster_column)
  batch_sym <- rlang::sym(batch_column)
  
  data <- seurat_obj@meta.data %>%
    dplyr::select(!!sym(cluster_column), !!sym(batch_column)) %>%
    dplyr::group_by(!!sym(cluster_column), !!sym(batch_column)) %>%
    dplyr::summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!sym(batch_column), values_from = count, values_fill = 0)
  
  if (normalize){
    
    batch_names <- setdiff(names(data), cluster_column)
    
    data_transformed <- data %>%
      dplyr::mutate(across(all_of(batch_names), ~ . / mean(.))) %>%
      dplyr::mutate(!!cluster_sym := data[[cluster_column]])
    
  } else{
    
    data_transformed <- data
    
  }
  
  # Calculate the fractions to ensure each cluster sums to 1 (100%)
  batch_frac <- data_transformed %>%
    rowwise() %>%
    mutate(rowSum = sum(c_across(-!!sym(cluster_column)))) %>%
    ungroup() %>%
    mutate(across(-c(!!sym(cluster_column)), ~ . / rowSum)) %>%
    dplyr::select(!rowSum)
  
  # Reshape for ggplot
  plot_data <- batch_frac %>%
    tidyr::pivot_longer(cols = -!!sym(cluster_column), names_to = batch_column, values_to = "fraction") 

  # Plot the stacked bar chart
  ggplot(plot_data, aes(x = factor(!!sym(cluster_column)), y = fraction, fill = !!sym(batch_column))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = cluster_column, y = "Fraction", fill = batch_column) +
    scale_y_continuous() + # Show percentages on y-axis - labels = scales::percent
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = palette)
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

#-----------------------------------------------------------------PATHS------------------------------------------------------

figures_folder <- '/projects/sle_jul_23_gabibov/luad/aucell/figures/'

#-----------------------------------------------------------------COLORS------------------------------------------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))


#-----------------------------------------------------------------READ DATA------------------------------------------------------

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated.rds')

integrated@meta.data <- integrated@meta.data %>%
  mutate(orig.ident_custom = case_when(
    stringr::str_starts(orig.ident, 'P') ~ 'WetLab data',
    TRUE ~ orig.ident
  ))

th1_geneset <- read.csv('/projects/sle_jul_23_gabibov/luad/aucell/genesets/GSE14308_TH1_VS_NAIVE_CD4_TCELL_UP.v2024.1.Hs.tsv', sep = '\t') %>%
  filter(STANDARD_NAME == 'GENE_SYMBOLS') %>% pull('GSE14308_TH1_VS_NAIVE_CD4_TCELL_UP')

th1_geneset <- unique(c("EOMES", "RUNX1", "RARA", "GZMA", "GZMH", "GZMB", "CCL4", "CCL5", "IFNG", "CXCR3", "BHLHE40", "GZMK", "EOMES", 
                 "RUNX1", "RARA", "GZMA", "GZMH", "GZMB", "CCL4", "CCL5", "IFNG", "CXCR3", "NKG7", "PRF1", "GNLY", "GZMA", 
                 "GZMH", "GZMB", "IFNG", "CXCR3", "BHLHE40", "KLRG1", "KLRD1", "KLRF1", "CX3CR1", "S1PR5", "CCL3", "CCL4", "CCL5"))

#-----------------------------------------------------------------PREPARE RAW DATA------------------------------------------------------

integrated_kras <- subset(
  integrated,
  subset = KRAS != "Unknown" & CelltypeAnnotationHigh == "Tem Th1 like" & (is.na(Treatment) | Treatment == 'pre'))

integrated_kras_mut <- GetAssayData(object = subset(
  integrated,
  subset = KRAS == "mut" & CelltypeAnnotationHigh == "Tem Th1 like" & (is.na(Treatment) | Treatment == 'pre')),
  layer = 'counts')


integrated_kras_wt <- GetAssayData(object = subset(
  integrated,
  subset = KRAS == "WT" & CelltypeAnnotationHigh == "Tem Th1 like" & (is.na(Treatment) | Treatment == 'pre')),
  layer = 'counts')

integrated_kras_mut <- GetAssayData(object = subset(
  integrated,
  subset = KRAS == "mut" & CelltypeAnnotationHigh == "Tem Th1 like" & (is.na(Treatment) | Treatment == 'pre')),
  layer = 'counts')

#-----------------------------------------------------------------CALCULATE CELL RANKINGS------------------------------------------------------

cell_ranrikgs_integrated_kras_wt <- AUCell::AUCell_buildRankings(integrated_kras_wt)

cell_ranrikgs_integrated_kras_mut <- AUCell::AUCell_buildRankings(integrated_kras_mut)


#-----------------------------------------------------------------CALCULATE CELL RANKINGS------------------------------------------------------

integrated_kras <- AddModuleScore(
  object = integrated_kras,
  features = list(th1_geneset),
  name = "Th1_geneset"
)

integrated_kras_meta <- integrated_kras@meta.data

ggplot(data = integrated_kras_meta, aes(x = KRAS, y = Th1_geneset1, fill = KRAS)) +
  geom_boxplot(outliers = FALSE) +
  #geom_point(position = position_dodge(width = 0.75), size = 1, aes(fill = KRAS)) +
  xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        panel.grid.major.x = element_blank())



celltypes_mapping_thesis <- c(
  '0_0' = "TcmI",
  '0_1' = "Trm",
  '0_2' = "Tem Th1 like GZMK",
  '0_3' = "Tcm Th1 like", # TGF-beta high
  '0_4' = "Tem Th1 like",
  '1' = "Tem T reg",
  '2_0' = "Tem Th1 like", 
  '2_1' = "Tem Tfh like", 
  '2_2' = "Exhausted PD1High Th1 like",
  '2_3' = "Exhausted PD1High Th1 like", 
  '2_4' = "Tem Tfh like",
  '2_5' = "Exhausted PD1High Th1 like", 
  '3' = "Tem T reg",
  '4_0' = "Tem Th1-Th17 like",
  '4_1' = "TcmI",
  '4_2' = "Tem Tfh like",
  '4_3' = "Tem Th1 like", # сомнительно
  '4_4' = "Tem Th1-Th17 like",
  '4_5' = "IFN-induced",
  '5_0' = "Unknown",
  '5_1' = "Trm",
  '5_2' = "TcmI",
  '5_3' = "Unknown",
  '5_4' = "TcmI",
  '5_5' = "Unknown",
  '6_0' = "Unknown",
  '6_1' = "TcmI",
  '6_2' = "Unknown",
  '6_3' = "Tem Th1 like", # High TNF level
  '6_4' = "Tem Th1 like",
  '6_5' = "Tem Th1 like", # High TNF level
  '7' = "TcmI",
  '8_0' = "Tem Th1 like",
  '8_1' = "Unknown",
  '8_2' = "TcmI", 
  '9' = "Tem Th1 like",
  '10' = "Tem Tfh like",
  '11_0' = "Tem Th17 like",
  '11_1' = "Tem Th17 like",
  '11_2' = "Unknown",
  '11_3' = "Unknown",
  '12' = "Tem Th1 like GZMK",
  '13' = "Tem T reg",
  '14_0' = "Tem T reg",
  '14_1' = "Tem Th1 like GZMK",
  '14_2' = "Tem Tfh like",
  '14_3' = "Tem Th1 like GZMK",
  '14_4' = "Unknown",
  '14_5' = "Tem T reg",
  '14_6' = "Tem Tfh like",
  '14_7' = "Tem T reg",
  '15' = "Tem Th17 like",
  '16' = "Tem T reg", # Maybe Tcm T reg??
  '17_0' = "Tem T reg",
  '17_1' = "Tem T reg",
  '17_2' = "Tem T reg",
  '17_3' = "Tem T reg",
  '18_0' = "Unknown",
  '18_1' = "TcmII",
  '18_2' = "Tem T reg",
  '18_3' = "Unknown",
  '18_4' = "Tcm Th17 like",
  '18_5' = "TcmII",
  '18_6' = "Tem Tfh like",
  '18_7' = "Tem Th1 like",
  '19' = "Tem T reg",
  '20' = "IFN-induced",
  '21_0' = "Tem Tfh like",
  '21_1' = "PD1High",
  '21_2' = "PD1High",
  '21_3' = "Tem Tfh like",
  '22' = "Temra",
  '23_0' = "Unknown",
  '23_1' = "Unknown",
  '23_2' = "Unknown",
  '23_3' = "Tem Th1-Th17 like",
  '23_4' = "Unknown",
  '23_5' = "Trm",
  '24' = "TcmI",
  '25' = "Cycling",
  '26' = "Cycling",
  '27_0' = "Unknown",
  '27_1' = "Unknown", # High TNF level
  '28_0' = "Unknown",
  '28_1' = "Unknown",
  '29' = "Unknown"
)

integrated@meta.data$CelltypeAnnotationThesis <- celltypes_mapping_thesis[integrated@meta.data$subclusters_final]

integrated@meta.data$CelltypeAnnotationThesis <- factor(integrated@meta.data$CelltypeAnnotationThesis, levels = c('TcmI',
                                                                                                                  'TcmII',
                                                                                                                  'Trm',
                                                                                                                  'Tcm Th17 like',
                                                                                                                  'Tem Th17 like',
                                                                                                                  'Tem Th1-Th17 like',
                                                                                                                  'Tcm Th1 like',
                                                                                                                  'Tem Th1 like',
                                                                                                                  'Tem Th1 like GZMK',
                                                                                                                  'Temra',
                                                                                                                  'Tem Tfh like',
                                                                                                                  'PD1High',
                                                                                                                  'Exhausted PD1High Th1 like',
                                                                                                                  'Tem T reg',
                                                                                                                  'IFN-induced',
                                                                                                                  'Cycling',
                                                                                                                  'Unknown'))



gse162500 <- readRDS('/projects/sle_jul_23_gabibov/luad/integrated_rds/gse162500/integrated_cd4_t_cell_data_vdj.rds')
gse162500@assays$RNA@layers$scale.data <- NULL

abnormal_cells_to_exclude_gse243013 <- readLines('/projects/sle_jul_23_gabibov/luad/integration_big/abnormal_cell_to_exclude_gse243013.txt')

cells_to_include <- setdiff(Cells(gse243013), abnormal_cells_to_exclude_gse243013)

gse243013 <- readRDS('/projects/sle_jul_23_gabibov/luad/integrated_rds/gse243013/integrated_cd4_t_cell_data_vdj.rds')
#gse243013@assays$RNA@layers$scale.data <- NULL

gse243013 <- subset(gse243013, cells = cells_to_include)

data_new <- merge(x = gse162500, y = gse243013)







