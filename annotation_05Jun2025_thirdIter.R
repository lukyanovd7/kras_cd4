.libPaths(c( .libPaths(), "/home/ssuleimanov/R/x86_64-pc-linux-gnu-library/4.1")) 

#-----------------------------------------------------------------LIBRARIES------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
library(cowplot)
library(patchwork)

setwd('/projects/sle_jul_23_gabibov/')

#-----------------------------------------------------------------COLORS------------------------------------------------------

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'

#-----------------------------------------------------------------READ DATA------------------------------------------------------

integration_methods <- 'rpca'

path_to_data <- paste0('/projects/sle_jul_23_gabibov/luad/integration_big_thirdIter/', integration_methods, '/luad_', integration_methods, '.rds')

integrated <- readRDS(path_to_data)

figures_folder <- "/projects/sle_jul_23_gabibov/luad/annotation_big/figures/05Jun25_thirdIter/"
tables_folder <- "/projects/sle_jul_23_gabibov/luad/annotation_big/tables/05Jun25_thirdIter/"

#-----------------------------------------------------------------ANNOTATE WITH PREVIOUS ANNOTATION (FIRST VERSION)------------------------------------------------------

annotation_df <- read.csv('/projects/sle_jul_23_gabibov/luad/annotation_big/data/luad_rpca_annotated_17Apr25.csv', row.names = 1)

annotation_df <- annotation_df[colnames(integrated), ]

integrated <- AddMetaData(integrated, select(annotation_df, CelltypeAnnotationLow))

integrated@meta.data$CelltypeAnnotationLow <- factor(integrated_ann@meta.data$CelltypeAnnotationLow, levels = names(subsets_cols))

DimPlot(integrated_ann, group.by='CelltypeAnnotationLow', cols=subsets_cols, raster = FALSE, label = FALSE)
ggsave2("previous_annotation_umap.png", path = figures_folder, width = 30, height = 15, units = "cm")

#-----------------------------------------------------------------PLOT CLUSTER 4------------------------------------------------------

integrated@meta.data$cluster4 <- ifelse(integrated@meta.data$seurat_clusters == 4, TRUE, FALSE)

DimPlot(integrated, group.by='cluster4', cols=c('lightgrey', 'red'), raster = FALSE, label = FALSE)
ggsave2("umap_cluster4.png", path = figures_folder, width = 20, height = 15, units = "cm")

#-----------------------------------------------------------------ANNOTATE WITH PREVIOUS ANNOTATION (FIRST VERSION)------------------------------------------------------

integrated_meta <- integrated@meta.data

integrated_subset_distribution <- integrated_meta %>% 
  filter(!is.na(CelltypeAnnotationLow)) %>% 
  select(c(seurat_clusters, CelltypeAnnotationLow)) %>% 
  group_by(seurat_clusters) %>% 
  count(CelltypeAnnotationLow, name='Count') %>% 
  mutate(MajorityVotingPrediction = Count / sum(Count)*100) %>% 
  arrange(seurat_clusters, desc(MajorityVotingPrediction)) %>% 
  ungroup()

integrated_subset_distribution %>% 
  group_by(seurat_clusters) %>% 
  mutate(Top = dense_rank(desc(MajorityVotingPrediction))) %>% 
  ungroup() %>% 
  mutate(Top=as.factor(Top)) %>%
  ggplot(aes(y=Top, x=MajorityVotingPrediction, fill=Top)) +
  geom_boxplot() +
  geom_point(size=2, position='jitter') +
  theme_bw() +
  xlim(0,100) +
  ggtitle('Subset frequency in annotated data') +
  scale_fill_manual(values=col_vector)

ggsave2("majority_voting.png", path = figures_folder, width = 15, height = 15, units = "cm")

integrated_top1_subset <- integrated_meta %>% 
  filter(!is.na(CelltypeAnnotationLow)) %>% 
  select(c(seurat_clusters, CelltypeAnnotationLow)) %>% 
  group_by(seurat_clusters) %>% 
  count(CelltypeAnnotationLow, name='Count') %>% 
  mutate(MajorityVotingPrediction = Count / sum(Count)*100) %>% 
  slice_max(MajorityVotingPrediction, n = 1, with_ties = FALSE) %>% 
  ungroup()

#-----------------------------------------------------------------FIND SUBCLUSTERS------------------------------------------------------

subclustering <- c(2, 4, 16, 18)

for (cluster in subclustering){
  integrated <- FindSubCluster(integrated, cluster = cluster, subcluster.name = paste0(cluster, '_sub'), graph.name = 'RNA_snn', resolution = 0.25)
}

integrated_meta <- integrated@meta.data

integrated_meta_subclusters <- select(integrated_meta, contains('_sub'))

result <- apply(integrated_meta_subclusters, 1, function(row) {
  unique_vals <- unique(row)
  if (length(unique_vals) == 1) {
    return(unique_vals)
  } else {
    tab <- table(row)
    unique_val <- names(tab[tab == 1])
    return(unique_val)
  }
})

integrated_meta$subclusters_final <- result
integrated@meta.data$subclusters_final <- result

#-----------------------------------------------------------------ANNOTATE SUBCLUSTERS WITH PREVIOUS ANNOTATION (FIRST VERSION)------------------------------------------------------

integrated_meta <- integrated@meta.data

integrated_subset_distribution <- integrated_meta %>% 
  filter(!is.na(CelltypeAnnotationLow)) %>% 
  select(c(subclusters_final, CelltypeAnnotationLow)) %>% 
  group_by(subclusters_final) %>% 
  count(CelltypeAnnotationLow, name='Count') %>% 
  mutate(MajorityVotingPrediction = Count / sum(Count)*100) %>% 
  arrange(subclusters_final, desc(MajorityVotingPrediction)) %>% 
  ungroup()

# integrated_subset_distribution %>% 
#   group_by(subclusters_final) %>% 
#   mutate(Top = dense_rank(desc(MajorityVotingPrediction))) %>% 
#   ungroup() %>% 
#   mutate(Top=as.factor(Top)) %>%
#   ggplot(aes(y=Top, x=MajorityVotingPrediction, fill=Top)) +
#   geom_boxplot() +
#   geom_point(size=2, position='jitter') +
#   theme_bw() +
#   xlim(0,100) +
#   ggtitle('Subset frequency in annotated data') +
#   scale_fill_manual(values=col_vector)
# 
# ggsave2("majority_voting.png", path = figures_folder, width = 15, height = 15, units = "cm")

integrated_top2_subset <- integrated_meta %>% 
  filter(!is.na(CelltypeAnnotationLow)) %>% 
  select(c(subclusters_final, CelltypeAnnotationLow)) %>% 
  group_by(subclusters_final) %>% 
  count(CelltypeAnnotationLow, name='Count') %>% 
  mutate(MajorityVotingPrediction = Count / sum(Count)*100) %>% 
  slice_max(MajorityVotingPrediction, n = 2, with_ties = FALSE) %>% 
  ungroup()

#-----------------------------------------------------------------EXPRESSION OF MARKER GENES------------------------------------------------------

gene_signatures <- list(
  "Tcm_I" = c("SELL", "TCF7", "IL7R", "CCR7"),
  "Tcm_II" = c("IL7R", "CCR7", "ANXA1", "ANXA2", "LMNA", "LGALS3", "LGALS1"),
  "Trm" = c("CD69", "NR4A1", "NR4A2", "ITGAE", "ITGA1", "DUSP6", "S1PR1", "ZNF683", 'PRDM1'), # https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.616309/full
  "Tem_Th17" = c("IL17A", "IL17F", "IL22", "IL26", "CCR6", "TNF", "IL23R", "RORC"),
  "Tem_Th1" = c("EOMES", "TBX22", "RUNX1", "RARA", "GZMA", "GZMH", "GZMB", "CCL4", "CCL5", "IFNG", "CXCR3", "CXCR4", "BHLHE40"),
  "Tem_Th2a" = c("GATA3", "STAT6", "IL4", "IL5", "IL13", "PTGDR2"), 
  "Tem_Th22" = c("CCR10"),
  "Tem_Th1_GZMK" = c("GZMK", "EOMES", "RUNX1", "RARA", "GZMA", "GZMH", "GZMB", "CCL4", "CCL5", "IFNG", "CXCR3", "CXCR4"),
  "Temra_Th1" = c("NKG7", "PRF1", "GNLY", "GZMA", "GZMH", "GZMB", "IFNG", "CXCR3", "BHLHE40", "KLRG1", "KLRD1", "KLRF1", "CX3CR1", "S1PR5", "CCL3", "CCL4", "CCL5"),
  "Tfh" = c("CXCR5", "BCL6", "ICA1", "IL6ST", "MAGEH1", "BTLA", "CD200", "CXCR3", "IL21", "ICOS", "PDCD1"),
  "Exhausted" = c("PDCD1", "LAG3", "TOX", "TOX2", "ICOS", "CXCL13", "TIGIT", "CTLA4", "HAVCR2"),
  "Treg_CD25low" = c("FOXP3", "FOXO1", "FOXO3", "CTLA4", "IL2RA", "IKZF2"),
  "Treg_CD25high" = c("FOXP3", "FOXO1", "FOXO3", "CTLA4", "IL2RA", "IKZF2"),
  "Tem_IFN-induced" = c("IFIT1", "IFIT3", "IFIT2", "STAT1", "MX1", "IRF7"),
  "Tcycle" = c("CDK1", "MKI67", "TUBA1B", "STMN1"),
  "Additional" = c('TGFB1')
)

gene_signatures_short <- c("MKI67", "TOP2A", "BIRC5", "IL22", "AHR", "CCR6", "CCR10", "PDCD1", "LAG3", "TIGIT", 
                           "CTLA4", "TOX", "BCL6", "CXCR5", "PDCD1", "ICOS", "IL21", "FOXP3", "IL2RA", "CTLA4", 
                           "TIGIT", "IKZF2", "TBX21", "IFNG", "GZMA", "GZMB", "PRF1", "TBX21", "IFNG", "IL12RB2", 
                           "CXCR3", "TCF7", "LEF1", "SELL", "CCR7", "IL7R", "ISG15", "IFIT1", "IFIT3", "MX1", 
                           "OAS1", "CCR7", "TCF7", "LEF1", "CD27", "IL7R", "GATA3", "IL4", "IL5", "IL13", "TCF7", 
                           "LEF1", "SELL", "CCR7", "IL7R", "PTCRA", "RORC", "IL17A", "IL17F", "IL22", "CCR6", 
                           "IL23R", "GATA3", "IL4", "IL5", "IL13", "TCF7", "LEF1", "SELL", "CCR7", "IL7R", "ISG15", "IFIT1", "MX1")

genes <- unique(unlist(gene_signatures))
cell_types <- rep(names(gene_signatures), sapply(gene_signatures, length))
gene_df <- data.frame(Cell_Type = cell_types, Gene = unlist(gene_signatures))

DotPlot(integrated, features = unique(gene_df$Gene), group.by = 'subclusters_final') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave2("marker_genes_subclustered.png", path = figures_folder, width = 45, height = 25, units = "cm")
ggsave2("marker_genes_subclustered.pdf", path = figures_folder, width = 45, height = 25, units = "cm")

# Subclusters for cluster 4

integrated@meta.data$subclusters_final_UMAP <- ifelse(integrated@meta.data$seurat_clusters %in% c(4), integrated@meta.data$subclusters_final, FALSE)

colors_umap_subclusters <- col_vector[1:10]
colors_umap_subclusters[4] <- 'lightgrey'

DimPlot(integrated, group.by = 'subclusters_final_UMAP', cols = colors_umap_subclusters, raster = FALSE)

ggsave2("subclusters_4_umap.png", path = figures_folder, width = 15, height = 15, units = "cm")

# Subclusters for cluster 2

integrated@meta.data$subclusters_final_UMAP <- ifelse(integrated@meta.data$seurat_clusters %in% c(2), integrated@meta.data$subclusters_final, FALSE)

colors_umap_subclusters <- col_vector[1:10]
colors_umap_subclusters[6] <- 'lightgrey'

DimPlot(integrated, group.by = 'subclusters_final_UMAP', cols = colors_umap_subclusters, raster = FALSE)

ggsave2("subclusters_2_umap.png", path = figures_folder, width = 15, height = 15, units = "cm")

#-----------------------------------------------------------------EXPRESSION OF MARKER GENES IN UMAP------------------------------------------------------

umap_features <- c('CCL5', 'GZMA', 'IFNG', 'IL17A', 'IL23R', 'RORC', 'CXCR5', 'IL6ST', 'BCL6', 'CXCL13', 'PDCD1', 'LAG3', 'CTLA4', 'TIGIT', 'TOX')

p_cluster <- DimPlot(integrated, group.by = "seurat_clusters", raster = FALSE, cols = col_vector, label = TRUE)

feature_plots <- lapply(umap_features, function(feature) {
  FeaturePlot(integrated, features = feature, raster = FALSE)
})

all_plots <- wrap_plots(c(list(p_cluster), feature_plots), ncol = 4)

ggsave2(filename = "combined_umap_plots.png",
        path = figures_folder,
        plot = all_plots,
        width = 25,
        height = 26,
        dpi = 300)


#-----------------------------------------------------------------ANNOTATION------------------------------------------------------

integrated@meta.data$CelltypeAnnotationLow <- NULL

celltypes_mapping_level1 <- c(
  '0' = 'T reg',
  '1' = 'T reg',
  '2_0' = 'PD1 high',
  '2_1' = 'Cytotoxic PD1 high',
  '2_2' = 'PD1 high',
  '2_3' = 'PD1 high',
  '2_4' = 'PD1 high',
  '3' = 'Tem Th1 like',
  '4_0' = 'Unknown',
  '4_1' = 'Tem Th17 like',
  '4_2' = 'Tem Th17 like',
  '5' = 'Tcm/Naive',
  '6' = 'T reg',
  '7' = 'Tem Th1 like',
  '8' = 'Trm Th1 like',
  '9' = 'T reg',
  '10' = 'Tem Tfh like',
  '11' = 'Tem Th1 like',
  '12' = 'Tem Th1 like',
  '13' = 'Tem Th1 like',
  '14' = 'Trm Th1 like',
  '15' = 'Tem Tfh like',
  '16_0' = 'Trm Tfh like',
  '16_1' = 'Trm Tfh like',
  '16_2' = 'Trm Tfh like',
  '17' = 'Cycling',
  '18_0' = 'Unknown',
  '18_1' = 'Unknown',
  '18_2' = 'Unknown',
  '19' = 'Cytotoxic TEMRA',
  '20' = 'IFN-induced',
  '21' = 'Tcm/Naive',
  '22' = 'T reg',
  '23' = 'Unknown'
)

celltypes_mapping_level2 <- c(
  '0' = 'T reg CD25 low',
  '1' = 'T reg CD25 high',
  '2_0' = 'PD1 high',
  '2_1' = 'Cytotoxic PD1 high',
  '2_2' = 'PD1 high',
  '2_3' = 'PD1 high',
  '2_4' = 'PD1 high',
  '3' = 'Tem Th1 like GZMA low',
  '4_0' = 'Unknown',
  '4_1' = 'Tem Th17 like CCR6 high IL17 low',
  '4_2' = 'Tem Th17 like CCR6 high IL17 high',
  '5' = 'Tcm/Naive',
  '6' = 'T reg CD25 low',
  '7' = 'Tem Th1 like GZMA low',
  '8' = 'Trm Th1 like TNF high',
  '9' = 'T reg CD25 low',
  '10' = 'Tem Tfh like CXCR5 low CXCL13 high',
  '11' = 'Tem Th1 like GZMA high',
  '12' = 'Tem Th1 like GZMK high',
  '13' = 'Tem Th1 like GZMA low',
  '14' = 'Trm Th1 like TNF high',
  '15' = 'Tem Tfh like CXCR5 high CXCL13 high',
  '16_0' = 'Trm Tfh like CXCR5 high CXCL13 low',
  '16_1' = 'Trm Tfh like CXCR5 high CXCL13 low',
  '16_2' = 'Trm Tfh like CXCR5 high CXCL13 low',
  '17' = 'Cycling',
  '18_0' = 'Unknown',
  '18_1' = 'Unknown',
  '18_2' = 'Unknown',
  '19' = 'Cytotoxic TEMRA',
  '20' = 'IFN-induced',
  '21' = 'Tcm/Naive',
  '22' = 'T reg CD25 low',
  '23' = 'Unknown'
)

celltypes_mapping_functional <- c(
  '0' = 'T reg',
  '1' = 'T reg',
  '2_0' = 'PD1 high',
  '2_1' = 'Cytotoxic PD1 high',
  '2_2' = 'PD1 high',
  '2_3' = 'PD1 high',
  '2_4' = 'PD1 high',
  '3' = 'Th1 like',
  '4_0' = 'Unknown',
  '4_1' = 'Th17 like',
  '4_2' = 'Th17 like ',
  '5' = 'Tcm/Naive',
  '6' = 'T reg',
  '7' = 'Th1 like',
  '8' = 'Th1 like',
  '9' = 'T reg',
  '10' = 'Tfh like',
  '11' = 'Th1 like',
  '12' = 'Th1 like',
  '13' = 'Th1 like',
  '14' = 'Th1 like',
  '15' = 'Tfh like',
  '16_0' = 'Tfh like',
  '16_1' = 'Tfh like',
  '16_2' = 'Tfh like',
  '17' = 'Cycling',
  '18_0' = 'Unknown',
  '18_1' = 'Unknown',
  '18_2' = 'Unknown',
  '19' = 'Cytotoxic TEMRA',
  '20' = 'IFN-induced',
  '21' = 'Tcm/Naive',
  '22' = 'T reg',
  '23' = 'Unknown'
)

integrated@meta.data$CelltypeAnnotationLevel1 <- celltypes_mapping_level1[integrated@meta.data$subclusters_final]

integrated@meta.data$CelltypeAnnotationLevel2 <- celltypes_mapping_level2[integrated@meta.data$subclusters_final]

integrated@meta.data$CelltypeAnnotationFunctional <- celltypes_mapping_functional[integrated@meta.data$subclusters_final]

integrated@meta.data$CelltypeAnnotationLevel1 <- factor(integrated@meta.data$CelltypeAnnotationLevel1, levels = c('Tcm/Naive',
                                                                                                            'Tem Th17 like',
                                                                                                            'Trm Th1 like',
                                                                                                            'Tem Th1 like',
                                                                                                            'Cytotoxic TEMRA',
                                                                                                            'Cytotoxic PD1 high',
                                                                                                            'Trm Tfh like',
                                                                                                            'Tem Tfh like',
                                                                                                            'PD1 high',
                                                                                                            'T reg',
                                                                                                            'IFN-induced',
                                                                                                            'Cycling',
                                                                                                            'Unknown'))

integrated@meta.data$CelltypeAnnotationLevel2 <- factor(integrated@meta.data$CelltypeAnnotationLevel2, levels = c('Tcm/Naive',
                                                                                                                  'Tem Th17 like CCR6 high IL17 low',
                                                                                                                  'Tem Th17 like CCR6 high IL17 high',
                                                                                                                  'Trm Th1 like TNF high',
                                                                                                                  'Tem Th1 like GZMA low',
                                                                                                                  'Tem Th1 like GZMA high',
                                                                                                                  'Tem Th1 like GZMK high',
                                                                                                                  'Cytotoxic TEMRA',
                                                                                                                  'Cytotoxic PD1 high',
                                                                                                                  'Trm Tfh like CXCR5 high CXCL13 low',
                                                                                                                  'Tem Tfh like CXCR5 low CXCL13 high',
                                                                                                                  'Tem Tfh like CXCR5 high CXCL13 high',
                                                                                                                  'PD1 high',
                                                                                                                  'T reg CD25 low',
                                                                                                                  'T reg CD25 high',
                                                                                                                  'IFN-induced',
                                                                                                                  'Cycling',
                                                                                                                  'Unknown'))

#-----------------------------------------------------------------REASSIGN KRAS STATUS FOR GSE243013------------------------------------------------------

integrated@meta.data$KRAS <- ifelse((integrated@meta.data$Driver.Mutations %in% c("RET", "BRAF", "ROS1", "HER2", "MET")) & (integrated@meta.data$orig.ident == 'GSE243013'), 'WT', integrated@meta.data$KRAS)


integrated@meta.data$Treatment <- ifelse(is.na(integrated@meta.data$Treatment), 'pre', integrated@meta.data$Treatment)
integrated@meta.data$Treatment <- ifelse(integrated@meta.data$Treatment == 'Post', 'post', integrated@meta.data$Treatment)

#-----------------------------------------------------------------MAKE UNIQUE PATIENT COLUMN------------------------------------------------------

integrated@meta.data$Patient <- paste0(integrated@meta.data$orig.ident, '_', integrated@meta.data$Patient)

length(unique(integrated@meta.data$Patient))

#-----------------------------------------------------------------INVERT UMAP2 COORDINATES-----------------------------------------------------

integrated[["umap"]]@cell.embeddings[, 2] <- integrated[["umap"]]@cell.embeddings[, 2] * -1

#-----------------------------------------------------------------HELPERS------------------------------------------------------

source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

#-----------------------------------------------------------------BASIC PLOTS AFTER ANNOTATION------------------------------------------------------

integrated@meta.data <- integrated@meta.data %>%
  mutate(orig.ident_custom = case_when(
    stringr::str_starts(orig.ident, 'P') ~ 'WetLab data',
    TRUE ~ orig.ident
  ))

integrated_clean <- subset(integrated, subset = CelltypeAnnotationLevel1 != 'Unknown')

DimPlot(integrated_clean, group.by='CelltypeAnnotationLevel1', cols=subsets_cols_level1, raster=FALSE) +
  ggtitle('CelltypeAnnotationLevel1') +
  xlab('UMAP1') +
  ylab('UMAP2')

ggsave2("clusters_umap_annotated_level1.png", path = figures_folder, width = 25, height = 19, units = "cm", dpi = 300)

DimPlot(integrated_clean, group.by='CelltypeAnnotationLevel2', cols=subsets_cols_level2, raster=FALSE) +
  ggtitle('CelltypeAnnotationLevel2') +
  xlab('UMAP1') +
  ylab('UMAP2')

ggsave2("clusters_umap_annotated_level2.png", path = figures_folder, width = 25, height = 19, units = "cm", dpi = 300)

plot_data <- plot_batch_fr(integrated_clean, cluster_column = 'CelltypeAnnotationLevel1', batch_column = 'orig.ident_custom', palette = col_vector, normalize = FALSE) +
  labs(fill = 'Datasets')
ggsave2("batch_disribution_over_annotation_level1.png", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)

plot_data <- plot_batch_fr(integrated_clean, cluster_column = 'CelltypeAnnotationLevel1', batch_column = 'orig.ident_custom', palette = col_vector, normalize = TRUE) +
  labs(fill = 'Datasets')
ggsave2("batch_disribution_over_annotation_level1_norm.png", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)

plot_data <- plot_batch_fr(integrated_clean, cluster_column = 'CelltypeAnnotationLevel2', batch_column = 'orig.ident_custom', palette = col_vector, normalize = FALSE) +
  labs(fill = 'Datasets')
ggsave2("batch_disribution_over_annotation_level2.png", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)

plot_data <- plot_batch_fr(integrated_clean, cluster_column = 'CelltypeAnnotationLevel2', batch_column = 'orig.ident_custom', palette = col_vector, normalize = TRUE) +
  labs(fill = 'Datasets')
ggsave2("batch_disribution_over_annotation_level2_norm.png", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)

#-----------------------------------------------------------------QUALITY CONTROL---------------------------------------------------------

quality_plot_custom(integrated_clean, group.by = 'CelltypeAnnotationLevel1', colors = subsets_cols_level1)

ggsave2("quality_control_by_annotation_level1.png", path = figures_folder, width = 20, height = 20, units = "cm", dpi = 600)

quality_plot_custom(integrated_clean, group.by = 'CelltypeAnnotationLevel2', colors = subsets_cols_level2)

ggsave2("quality_control_by_annotation_level2.png", path = figures_folder, width = 35, height = 25, units = "cm", dpi = 600)

quality_plot_custom(integrated_clean, group.by = 'orig.ident_custom', colors = col_vector)

ggsave2("quality_control_by_batch.png", path = figures_folder, width = 20, height = 20, units = "cm", dpi = 600)

#-----------------------------------------------------------------SUBSAMPLE EACH BATCH TO 5K CELLS------------------------------------------------------

subsampling_number <- 5000

cells_to_keep <- integrated_clean@meta.data %>%
  tibble::rownames_to_column(var = 'barcodes') %>%
  group_by(orig.ident_custom) %>%
  group_map(~ {
    cells_in_batch <- (.x$barcodes)
    n <- length(cells_in_batch)
    if (n > subsampling_number) {
      sample(cells_in_batch, subsampling_number)
    } else {
      cells_in_batch
    }
  }) %>%
  unlist()

integrated_clean_subsampled <- subset(integrated_clean, cells = cells_to_keep)

DimPlot(integrated_clean_subsampled, group.by='orig.ident_custom', cols=col_vector, raster=FALSE) +
  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position = 'bottom')

ggsave2("datasets_split_umap_subsampled5k.png", path = figures_folder, width = 30, height = 34, units = "cm", dpi = 600)


DimPlot(integrated_clean_subsampled, group.by='CelltypeAnnotationLevel1', split.by = 'orig.ident_custom', cols=subsets_cols_level1, raster=FALSE, pt.size = 1.5) +
  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position = 'bottom')

ggsave2("datasets_split_umap_subsampled5k_annotation_level1.png", path = figures_folder, width = 35, height = 40, units = "cm", dpi = 600)

DimPlot(integrated_clean_subsampled, group.by='CelltypeAnnotationLevel2', split.by = 'orig.ident_custom', cols=subsets_cols_level2, raster=FALSE, pt.size = 1.5) +
  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position = 'bottom')

ggsave2("datasets_split_umap_subsampled5k_annotation_level2.png", path = figures_folder, width = 35, height = 40, units = "cm", dpi = 600)


#-----------------------------------------------------------------SAVE ANNOTATED OBJECT------------------------------------------------------

saveRDS(integrated, '/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')

#-----------------------------------------------------------------3D UMAP------------------------------------------------------

integrated_3Dumap <- RunUMAP(integrated_clean, reduction = "integrated.rpca", dims = 1:30, n.components = 3L)

plot.data <- FetchData(object = integrated_3Dumap, vars = c("umap_1", "umap_2", "umap_3", "CelltypeAnnotationLevel1", "CelltypeAnnotationLevel2"), slot = 'data')

plot.data$color_level1 <- subsets_cols_level1[plot.data$CelltypeAnnotationLevel1]
plot.data$color_level2 <- subsets_cols_level2[plot.data$CelltypeAnnotationLevel2]

palette_named_level1 <- setNames(unique(plot.data$color_level1), unique(plot.data$CelltypeAnnotationLevel1))
palette_named_level2 <- setNames(unique(plot.data$color_level2), unique(plot.data$CelltypeAnnotationLevel2))

plotly_graph_level1 <- plotly::plot_ly(
            data = plot.data,
            x = ~umap_1,
            y = ~umap_2,
            z = ~umap_3,
            type = "scatter3d",
            mode = "markers",
            color = ~CelltypeAnnotationLevel1,
            colors = palette_named_level1,
            marker = list(size = 5, opacity = 0.5),
            hoverinfo = "text"
)

plotly_graph_level2 <- plotly::plot_ly(
  data = plot.data,
  x = ~umap_1,
  y = ~umap_2,
  z = ~umap_3,
  type = "scatter3d",
  mode = "markers",
  color = ~CelltypeAnnotationLevel2,
  colors = palette_named_level2,
  marker = list(size = 5, opacity = 0.5),
  hoverinfo = "text"
)

plotly_graph_level1
plotly_graph_level2
