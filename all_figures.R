### Generating Figures


# Figure 1B and 1C UMAP and batch distribution--------------------------------------------------------


library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(mltest)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(mltest)
library(pheatmap)

source('/path_to_folder/colors.R')
source('/path_to_folder/functions.R')

figures_folder <- "/path_to_folder/"
tables_folder <- "/path_to_folder/"

working_folder <- 'rpca'

integrated <- readRDS('/path_to_folder/integrated.rds')

integrated_backup <- integrated

integrated@meta.data <- integrated@meta.data %>%
  mutate(orig.ident_custom = case_when(
    stringr::str_starts(orig.ident, 'P') ~ 'WetLab data',
    TRUE ~ orig.ident
  ))

integrated_clean <- subset(integrated, subset = CelltypeAnnotationLevel1  != 'Unknown')
plot_data <- plot_batch_fr(integrated, cluster_column = 'seurat_clusters', batch_column = 'orig.ident_custom', palette = col_vector, normalize = FALSE) +
  labs(fill = 'Datasets')

plot_data <- plot_batch_fr(integrated, cluster_column = 'seurat_clusters', batch_column = 'orig.ident_custom', palette = col_vector, normalize = TRUE) +
  labs(fill = 'Datasets')

subsets_cols <- c(
  'Tcm/Naive' = '#1B9E77',                       # Teal green
  'Trm' = '#01665E',                        # Dark cyan
  'Tcm Th1-Th17 like' = '#1F78B4',          # Deep blue
  'Trm Th17 like CCR6 low' = '#A6CEE3',     # Light sky blue
  'Tem Th17 like CCR6 high' = '#6A3D9A',    # Dark purple
  'Tem Th17 like CCR6 low' = '#CAB2D6',     # Light lavender
  'Trm Th1 like' = '#FF7F00',               # Bright orange
  'Tem Th1 like' = '#FDBF6F',               # Soft orange
  'Tem Th1 like GZMK' = '#B15928',          # Brownish-orange
  'Tem Th1 like Exhausted PD1 high' = '#E31A1C',  # Bright red
  'Temra' = '#FB9A99',                      # Light red/pink
  'Tcm Tfh like CXCL13 low' = '#B2DF8A',    # Light green
  'Trm Tfh like CXCL13 low' = '#33A02C',    # Medium green
  'Tem Tfh like CXCL13 high' = '#006400',   # Dark green
  'Tem Tfh like CXCL13 low' = '#8DD3C7',    # Aqua green
  'Exhausted PD1 high' = '#FDDA0D',         # Yellow
  'T reg CD25 high' = '#BEAED4',            # Light purple
  'T reg CD25 low' = '#984EA3',             # Dark purple
  'IFN-induced' = '#999999',                # Neutral gray
  'Cycling high' = '#E6AB02',               # Mustard yellow
  'Cycling low' = '#7D6231'                 # Dark brown
)

integrated@meta.data <- integrated@meta.data %>%
  mutate(orig.ident_custom = case_when(
    stringr::str_starts(orig.ident, 'P') ~ 'WetLab data',
    TRUE ~ orig.ident
  ))

integrated_clean <- subset(integrated, subset = CelltypeAnnotationLevel2 %in% setdiff(unique(integrated@meta.data$CelltypeAnnotationLevel2), c('Artefact', 'Unknown')))

DimPlot(integrated_clean, group.by='CelltypeAnnotationLevel2', cols=subsets_cols_level2, raster=FALSE) + 
  ggtitle('Integrated datasets') +
  xlab('UMAP1') +
  ylab('UMAP2')
#1B
ggsave2("1B_clusters_umap_annotated_new.pdf", path = figures_folder, width = 25, height = 19, units = "cm", dpi = 300)


#DimPlot(integrated_clean, group.by='orig.ident_custom', cols=col_vector, raster=FALSE) +
#  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
#  xlab('UMAP1') +
#  ylab('UMAP2') +
#  theme(legend.position = 'bottom')
#ggsave2("datasets_split_umap.png", path = figures_folder, width = 30, height = 30, units = "cm", dpi = 300)

plot_data <- plot_batch_fr(integrated_clean, cluster_column = 'CelltypeAnnotationLevel2', batch_column = 'orig.ident_custom', palette = col_vector, normalize = TRUE) +
  labs(fill = 'Datasets')
#1C
ggsave2("1C_normalized_batch_disribution_over_clusters_annotated_new.pdf", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)



# S3 Quality control----------------------------------------------------------------------

#S3
quality_plot_custom(integrated_clean, group.by = 'CelltypeAnnotationLevel2') 
ggsave2("S3_quality_control_by_clusters_annotated.pdf", path = figures_folder, width = 30, height = 34, units = "cm", dpi = 300)
quality_plot_custom(integrated_clean, group.by = 'orig.ident_custom', colors = col_vector) 
ggsave2("S1_quality_control_by_batches_annotated.pdf", path = figures_folder, width = 30, height = 34, units = "cm", dpi = 300)



# Potential figure Number of cells per patient -------------------------------------------------------


cell_counts <- integrated@meta.data %>%
  count(Patient, name = "CellCount")

cell_counts_kras <- integrated@meta.data %>%
  filter(KRAS != 'Unknown') %>%
  count(Patient, name = "CellCount") %>%
  left_join(integrated@meta.data %>% filter(KRAS != 'Unknown') %>% select(Patient, KRAS) %>% distinct(), by = "Patient")

ggplot(cell_counts, aes(x = CellCount)) +
  geom_histogram(color = "black", fill = "steelblue") +
  labs(title = "Distribution of Cell Counts per Patient",
       x = "Number of Cells",
       y = "Number of Patients") +
  theme_minimal()

ggplot(cell_counts_kras, aes(x = CellCount, fill = KRAS)) +
  geom_density(color = "black") +
  labs(title = "Distribution of Cell Counts per Patient",
       x = "Number of Cells",
       y = "Number of Patients") +
  theme_minimal()

ggplot(cell_counts_kras, aes(x = KRAS, y = CellCount, fill = KRAS)) +
  geom_boxplot(color = "black") +
  geom_point(position = 'jitter') +
  labs(title = "Distribution of Cell Counts per Patient",
       x = "Number of Cells",
       y = "Number of Patients") +
  theme_minimal()


# 1D Heatmap----------------------------------------------------------------------


Idents(integrated_clean) <- "CelltypeAnnotationLevel2"

gene_signatures <- list(
  "Tcm/Naive" = c("SELL", "TCF7", "CCR7", "LMNA"),
  "Trm" = c("NR4A1", "NR4A2"),
  "Th17" = c("IL17A", "IL17F", "IL22", "IL26", "CCR6", "IL23R", "RORC"),
  "Th1" = c("EOMES", "RUNX1", "RARA", "GZMA", "GZMH", "GZMB", "CCL3", "CCL4", "CCL5", "IFNG", "BHLHE40", "GZMK"), # "TBX22",
  "Temra" = c("NKG7", "PRF1", "GNLY", "KLRG1", "KLRD1", "KLRF1", "CX3CR1", "S1PR5"),
  "Tfh" = c("CXCR5", "BCL6", "ICA1", "IL6ST", "BTLA", "CD200", "IL21", "ICOS"),
  "Exhausted" = c("PDCD1", "LAG3", "TOX", "TOX2", "CXCL13", "TIGIT", "CTLA4", "HAVCR2"),
  "Treg" = c("FOXP3", "IL2RA", "IKZF2"),
  "IFN-induced" = c("IFIT1", "IFIT3", "IFIT2", "STAT1", "MX1", "IRF7"),
  "Tcycle" = c("CDK1", "MKI67", "STMN1")
)

gene_group_df <- do.call(rbind, lapply(names(gene_signatures), function(group) {
  data.frame(Gene = gene_signatures[[group]], Group = group)
}))

# Get average expression matrix
avg_exp <- AverageExpression(integrated_clean, assays = 'RNA', features = unique(gene_group_df$Gene), group.by = "CelltypeAnnotationLevel2", layer = "data")$RNA

# Define ordered cell types and their groupings
celltype_order <- c(
  "Tcm/Naive",                                               # Group 1
  "Tem Th17 like CCR6 high IL17 low",                        # Group 2
  "Tem Th17 like CCR6 high IL17 high",
  "Tem Th1 like GZMA low",                                   # Group 3
  "Tem Th1 like GZMA high",
  "Tem Th1 like GZMK high",
  "Trm Th1 like TNF high",
  "Cytotoxic TEMRA",                                         # Group 4
  "Cytotoxic PD1 high",                                      # Group 5
  "PD1 high",
  "Trm Tfh like CXCR5 high CXCL13 low",                      # Group 6
  "Tem Tfh like CXCR5 low CXCL13 high",
  "Tem Tfh like CXCR5 high CXCL13 high",
  "T reg CD25 low",                                          # Group 7
  "T reg CD25 high",
  "IFN-induced",                                             # Group 8
  "Cycling"                                                  # Group 9
)

gene_order <- rev(intersect(unique(gene_group_df$Gene), rownames(avg_exp)))

# Subset and reorder matrix
ordered_mat <- avg_exp[gene_order, celltype_order]

# Scale per gene (optional, improves contrast)
ordered_mat_scaled <- scale(t(ordered_mat))

# ANNOTATION CELLS
annotation_row <- data.frame(CellType = factor(celltype_order, levels = celltype_order))
rownames(annotation_row) <- celltype_order

annotation_row_colors <- list(CellType = subsets_cols_level2)

# ANNOTATION GENES

rownames(gene_group_df) <- gene_group_df$Gene

annotation_col <- data.frame(
  Phenotype = gene_group_df[colnames(ordered_mat_scaled), "Group"]
)

annotation_col$Phenotype <- factor(annotation_col$Phenotype, levels = names(gene_signatures))

rownames(annotation_col) <- colnames(ordered_mat_scaled)

group_levels <- unique(annotation_col$Phenotype)
annotation_colors_col <- list(
  Phenotype = setNames(
    brewer.pal(length(group_levels), "Set3")[seq_along(group_levels)],
    group_levels
  )
)

gaps_col <- annotation_col %>%
  dplyr::mutate(Row = row_number()) %>%
  dplyr::group_by(Phenotype) %>%
  dplyr::summarise(gap_position = max(Row)) %>%
  dplyr::pull(gap_position)

group_sizes <- c(1, 2, 4, 1, 2, 3, 2, 1, 1)

gaps_row <- cumsum(head(group_sizes, -1))

p1 <- pheatmap(as.matrix(ordered_mat_scaled),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               #annotation_row = annotation_row,
               annotation_col = annotation_col,
               annotation_colors = c(annotation_colors_col), # Add if needed - annotation_row_colors
               gaps_col = gaps_col,
               gaps_row = gaps_row,
               fontsize_row = 8,
               fontsize_col = 8,
               annotation_legend = TRUE,
               main = "",
               color = colorRampPalette(c("blue", "white", "#FDDA0D", "red"))(100),
               breaks = seq(-2, 4, length.out = 101),
               border_color = 'white')


pdf(file=paste0(figures_folder, "/1D_heatmap_nonclustered_horizontal.pdf"), width = 30/2.54, height = 9/2.54)
grid::grid.newpage()
grid::grid.draw(p1$gtable)
#draw(p1)
dev.off()

# Potential figure CLUSTER BY KRAS STATUS------------------------------------------------------

integrated_clean_kras <- subset(integrated_clean, subset = KRAS != 'Unknown')

plot_data <- plot_batch_fr(integrated_clean_kras, cluster_column = 'CelltypeAnnotationLevel2', batch_column = 'KRAS', palette = col_vector[c(19,30)], normalize = TRUE) +
  labs(fill = 'KRAS status')

ggsave2("norm_kras_disribution_over_celltypes.png", path = figures_folder, width = 20, height = 12, units = "cm", dpi = 300)


# S2 SUBSAMPLE EACH BATCH TO 5K CELLS ----------------------------------------------------------------------


# BATCHES BY KRAS STATUS
meta <- integrated_clean@meta.data
meta_kras <- integrated_clean_kras@meta.data

meta_kras %>% group_by(orig.ident_custom, KRAS) %>% summarize(n = n_distinct(Patient)) 

View(meta %>% group_by(orig.ident_custom, KRAS) %>% summarize(n = n_distinct(Patient)))

# SUBSAMPLE EACH BATCH TO 5K CELLS

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

#ggsave2("datasets_split_umap_subsampled5k.png", path = figures_folder, width = 30, height = 30, units = "cm", dpi = 300)

#S2
DimPlot(integrated_clean_subsampled, group.by='CelltypeAnnotationLevel2', split.by = 'orig.ident_custom', cols=subsets_cols_level2, raster=FALSE) +
  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position = 'bottom')

ggsave2("S2_datasets_split_umap_subsampled5k_celltypes.pdf", path = figures_folder, width = 30, height = 34, units = "cm", dpi = 300)

#S1 CELL_TYPES_QUALITY_METRIC------------------------------------------------------

#Idents(integrate) <- integrated@meta.data$CelltypeAnnotationLevel2

quality_plot_custom(integrated_clean, group.by = 'orig.ident_custom')
ggsave2("S1_CELL_TYPES_QUALITY_METRICS.pdf", path = figures_folder, width = 30, height = 34, units = "cm", dpi = 300)


# 1A Comparing integrations --------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(mltest)
library(ggpubr)
library(purrr)

base_path <- "/path_to_folder/integration_big_thirdIter"

methods <- c("rpca", "cca", "jpca", "harmony")

read_and_label <- function(method) {
  path <- file.path(base_path, method, "tables", paste0("majorityVoting_statistics_", method, ".csv"))
  read.csv(path) %>%
    mutate(method = method)
}

data <- map_dfr(methods, read_and_label)

median_labels <- data %>%
  group_by(method) %>%
  summarise(median_freq = median(MajorityVotingPrediction),
            xpos = 1)

median_labels_filtered <- data %>%
  filter((Frequency_in_SeuratCluster >= 1) & (Count >= 100)) %>%
  group_by(method) %>%
  summarise(median_freq = median(MajorityVotingPrediction),
            xpos = 1)

#-----------------------------------------------------------------SUBSET PREDICTION SCORE

p1 <- data %>% 
  ggplot(aes(y=method, x=MajorityVotingPrediction, fill=method)) +
  geom_boxplot() +
  geom_point(size=2, position='jitter') +
  xlim(0,100) +
  ggtitle('') +
  scale_fill_manual(values=col_vector) +
  ylab('Integration method') +
  xlab('Majority voting prediction') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid
    panel.grid.minor.y = element_blank(),  # Remove horizontal minor grid
    panel.grid.minor.x = element_blank()   # Optional: remove vertical minor grid
  ) +
  geom_text(data = median_labels,
            aes(x = xpos, y = method, label = round(median_freq, 2)),
            inherit.aes = FALSE,
            hjust = 0, size = 5)

p1

p2 <- data %>% 
  filter((Frequency_in_SeuratCluster >= 1) & (Count >= 100)) %>%
  ggplot(aes(y=method, x=MajorityVotingPrediction, fill=method)) +
  geom_boxplot() +
  geom_point(size=2, position='jitter') +
  xlim(0,100) +
  ggtitle('') +
  scale_fill_manual(values=col_vector) +
  ylab('Integration method') +
  xlab('Majority voting prediction') +
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_blank(),  # Remove horizontal major grid
    panel.grid.minor.y = element_blank(),  # Remove horizontal minor grid
    panel.grid.minor.x = element_blank()   # Optional: remove vertical minor grid
  ) +
  geom_text(data = median_labels_filtered,
            aes(x = xpos, y = method, label = round(median_freq, 2)),
            inherit.aes = FALSE,
            hjust = 0, size = 5)

p2
#1A
ggsave2("1A_comparison_integration_methods_without_abnormal.pdf", path = figures_folder, width = 12, height = 10, units = "cm", dpi=300)


# S4 heatmap top 5 --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)

source('/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/luad_figures_tables/colors_fixed.R')
source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')

subsets_cols <- subsets_cols_level2

#-----------------------------------------------------------------READ DATA-

figures_folder <- '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed'
data_folder <- '/projects/sle_jul_23_gabibov/luad/diff_expression/data/'

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')

#-----------------------------------------------------------------ASSIGN CORRECT VERSION OF ANNOTATION
integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel2

#-----------------------------------------------------------------DELETE UNKNOWN AND ARTEFACT CELL TYPES

integrated <- subset(integrated, subset = CelltypeAnnotationLevel1 != 'Unknown')

#-----------------------------------------------------------------HEATMAP

Idents(integrated) <- integrated@meta.data$CelltypeAnnotationLow

markers <- FindAllMarkers(integrated,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.58)

write.csv(markers, paste0(data_folder, 'degs_per_celltype.csv'))

# markers <- read.csv(paste0(data_folder, 'degs_per_celltype.csv'))

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

subsampling_number <- 500

cells_to_keep <- integrated@meta.data %>%
  tibble::rownames_to_column(var = 'barcodes') %>%
  group_by(CelltypeAnnotationLow) %>%
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

integrated_subsampled <- subset(integrated, cells = cells_to_keep)

Idents(integrated_subsampled) <- integrated_subsampled@meta.data$CelltypeAnnotationLow

DoHeatmap(integrated_subsampled, features = top5$gene, size = 4, group.colors = subsets_cols) +
  scale_fill_gradientn(colors = c("purple", "black", "yellow")) +
  theme(axis.text.x = element_blank()) + 
  guides(color = "none", fill = guide_colorbar(title = "Expression"))

cowplot::ggsave2('S4_heatmap5top.pdf', path = figures_folder, width = 30, height = 35, units = 'cm', dpi = 600)


# 2A Composition -------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))


figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras"
data_folder <- '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras'

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')


integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel2

#-----------------------------------------------------------------CELL TYPE COMPOSITION

PSTFX <- 'preposttreatment_mad_filtered_annotation_level2b'

celltype_composition <- integrated@meta.data %>%
  filter(KRAS != 'Unknown') %>%
  #filter(Treatment == 'pre') %>%
  #filter(Treatment == 'post') %>%
  filter(!CelltypeAnnotationLow %in% c('Unknown', 'Artefact')) %>%
  group_by(KRAS, Patient, CelltypeAnnotationLow) %>%
  summarise(Count = n()) %>%
  mutate(Total = sum(Count)) %>%
  filter(Total >= 200) %>%
  mutate(Frequency = (Count / Total)) %>%
  #filter(!CelltypeAnnotationLow %in% c('Unknown', 'Artefact')) %>%
  ungroup() %>%
  arrange(Patient, desc(Frequency))

# write.csv(celltype_composition, paste0(data_folder, 'composition_df_200cellF.csv'))

# CALCULATE MAD 

stats <- celltype_composition %>%
  group_by(CelltypeAnnotationLow, KRAS) %>%
  summarise(
    med = median(Frequency),
    mad_val = mad(Frequency),
    .groups = "drop"
  )

celltype_composition <- celltype_composition %>%
  left_join(stats, by = c("CelltypeAnnotationLow", "KRAS")) %>%
  mutate(
    mad_score = abs(Frequency - med) / mad_val
  )

# PLOT DISTRIBUTION OF MAD SCORES

threshold <- quantile(celltype_composition$mad_score, 0.9, na.rm = TRUE)

ggplot(celltype_composition, aes(x = mad_score)) +
  geom_histogram(bins = 105) +
  theme_bw() +
  xlab('MAD score') +
  ylab('Count') +
  geom_vline(xintercept = threshold, color = 'red', linetype = 'dashed')

ggsave2(paste0("madScore_", PSTFX, ".png"), path = figures_folder, width = 10, height = 8, units = "cm", dpi = 300)

# FILTER OUT ONLY SAMPLES WITH HIGH MAD

celltype_composition <- celltype_composition %>%
  filter(mad_score <= threshold)

# write.csv(celltype_composition, paste0(data_folder, 'composition_df_200cellF_madF.csv'))

####### STATISTICAL TESTING
wilcox_results <- celltype_composition %>%
  mutate(KRAS = as.factor(KRAS)) %>%
  group_by(CelltypeAnnotationLow) %>%
  summarise(
    p = tryCatch(
      wilcox.test(Frequency ~ KRAS, exact=FALSE)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  mutate(star = ifelse(p < 0.05, "", "")) %>%
  mutate(star_adj = ifelse(p_adj < 0.05, "*", "")) %>%
  mutate(star_sccoda = case_when(
    CelltypeAnnotationLow %in% c('Tem Th1 like GZMA high', 'Cytotoxic PD1 high', 'Tem Tfh like CXCR5 high CXCL13 high') ~ "",
    TRUE ~ ""
  ))

celltype_composition <- celltype_composition %>%
  left_join(wilcox_results, by = "CelltypeAnnotationLow")
####### STATISTICAL TESTING

# write.csv(wilcox_results, paste0(data_folder, 'celltype_composition_kras_preposttreatment_mad_filtered.csv'), row.names = FALSE)


patients_per_celltypes <- celltype_composition %>% group_by(CelltypeAnnotationLow, KRAS) %>% summarise(NumPatients = n_distinct(Patient), .groups = "drop")

y_limit <- ceiling(max(patients_per_celltypes$NumPatients) / 5) * 5

c1 <- ggplot(data = celltype_composition, aes(x = CelltypeAnnotationLow, y = Frequency, fill = KRAS))+
  geom_point(position = position_jitterdodge(), size = 1.5) +
  geom_boxplot(outliers = FALSE) +
  geom_text(
    data = wilcox_results,
    aes(y = max(celltype_composition$Frequency)+(max(celltype_composition$Frequency)/4), x = CelltypeAnnotationLow, label = star),
    inherit.aes = FALSE,
    vjust = 2,
    size = 5,
    color = 'red'
  ) +
  geom_text(
    data = wilcox_results,
    aes(y = max(celltype_composition$Frequency)+(max(celltype_composition$Frequency)/4), x = CelltypeAnnotationLow, label = star_adj),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 5,
    color = 'black'
  ) +
  geom_text(
    data = wilcox_results,
    aes(y = max(celltype_composition$Frequency)+(max(celltype_composition$Frequency)/4), x = CelltypeAnnotationLow, label = star_sccoda),
    inherit.aes = FALSE,
    vjust = 0.5,
    size = 5,
    color = 'navy'
  ) +
  xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 8),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = col_vector[c(19,30,40)])

c2 <- ggplot(data = celltype_composition, aes(x = CelltypeAnnotationLow, fill = KRAS)) +
  geom_bar(stat = 'count', position = 'dodge') +
  xlab('') +
  ylab('Patients') +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_y_continuous(
    limits = c(0, y_limit),
    breaks = seq(0, y_limit, by = 5)) +
  scale_fill_manual(values = col_vector[c(19,30,40)])

ggpubr::ggarrange(c2, c1, nrow = 2, align = 'v', common.legend = TRUE, legend = 'right', heights = c(1, 4), vjust = 0.1)

ggsave2(paste0("2A_celltype_composition_full_", PSTFX, ".pdf"), path = figures_folder, width = 24, height = 14, units = "cm", dpi = 600)

#-----------------------------------------------------------------SAVE FILTERED SEURAT

integrated@meta.data$Composition_Filtering <- paste0(integrated@meta.data$Patient, '_', integrated@meta.data$CelltypeAnnotationLow)

celltype_composition$Composition_Filtering <- paste0(celltype_composition$Patient, '_', celltype_composition$CelltypeAnnotationLow)

integrated_filtered <- subset(integrated, subset = Composition_Filtering %in% unique(celltype_composition$Composition_Filtering))

#saveRDS(integrated_filtered, paste0(data_folder, "celltype_composition_full_", PSTFX, ".rds"))

#write.csv(integrated_filtered@meta.data, paste0('/projects/sle_jul_23_gabibov/luad/seurat_to_anndata/metadata_full_', PSTFX,'.csv'))

rm(integrated_filtered, celltype_composition)

# Potential figure BARPLOTS WITH MUTATIONS------------------------------------------------------

integrated@meta.data %>%
  filter(KRAS != 'Unknown') %>%
  group_by(KRAS, Patient, Treatment) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(KRAS, Treatment) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(Treatment = paste0(Treatment, 'T')) %>%
  mutate(KRAS_Treatment = paste(KRAS, Treatment, sep = "-")) %>%
  ggplot(aes(x = KRAS, y = n, fill = KRAS_Treatment)) +
  geom_bar(stat = 'identity') +
  xlab('KRAS mutation status') +
  ylab('Count') +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "mut-preT" = "#B2DF8A",     # light green
      "mut-postT" = "#66A061",    # dark green
      "WT-preT" = "#B3CDE3",    # light blue
      "WT-postT" = "#4F81BD"    # dark blue
    )
  ) +
  labs(fill = 'Groups')

#ggsave2("number_of_patients.png", path = figures_folder, width = 8, height = 7, units = "cm", dpi = 300)

# Potential figure STABILITY OF COMPOSITIONAL SHIFT ACROSS FILTRATION BY AMOUNT OF CELLS---------

integrated <- subset(integrated, subset = CelltypeAnnotationLevel2 != 'Unknown')

number_of_cells <- c(200, 400, 600)

plot_list <- list()

for (celltype_of_interest in levels(integrated@meta.data$CelltypeAnnotationLevel2)) {
  benchmarking_list <- list()
  
  for (i in seq_along(number_of_cells)) {
    cells <- number_of_cells[i]
    
    celltype_composition <- integrated@meta.data %>%
      filter(KRAS != 'Unknown') %>%
      filter(!CelltypeAnnotationLow %in% c('Unknown', 'Artefact')) %>%
      group_by(KRAS, Patient, CelltypeAnnotationLow) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Patient) %>%
      mutate(Total = sum(Count)) %>%
      filter(Total > cells) %>%
      mutate(Frequency = Count / Total) %>%
      ungroup()
    
    stats <- celltype_composition %>%
      group_by(CelltypeAnnotationLow, KRAS) %>%
      summarise(
        med = median(Frequency),
        mad_val = mad(Frequency),
        .groups = "drop"
      )
    
    celltype_composition <- celltype_composition %>%
      left_join(stats, by = c("CelltypeAnnotationLow", "KRAS")) %>%
      mutate(mad_score = abs(Frequency - med) / mad_val) %>%
      filter(mad_score <= quantile(mad_score, 0.9, na.rm = TRUE)) %>%
      filter(CelltypeAnnotationLow == celltype_of_interest)
    
    benchmarking_list[[i]] <- celltype_composition
  }
  
  benchmarking_df <- bind_rows(
    lapply(seq_along(benchmarking_list), function(i) {
      df <- benchmarking_list[[i]]
      df$Filtering <- factor(paste0(number_of_cells[i], "_cells"))
      return(df)
    })
  )
  
  p_values <- benchmarking_df %>%
    group_by(Filtering) %>%
    summarise(
      p = tryCatch(
        wilcox.test(Frequency ~ KRAS, exact = FALSE)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      star = ifelse(p < 0.05, "*", ""),
      p_adj = p.adjust(p, method = "BH"),
      star_adj = ifelse(p_adj < 0.05, "*", "")
    )
  
  benchmarking_df <- benchmarking_df %>%
    left_join(p_values, by = "Filtering")
  
  c1 <- ggplot(benchmarking_df, aes(y = Filtering, x = Frequency, fill = KRAS)) +
    geom_point(position = position_jitterdodge(), size = 1.5) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(data = p_values,
              aes(x = max(benchmarking_df$Frequency, na.rm = TRUE) * 1.05,
                  y = Filtering,
                  label = star),
              inherit.aes = FALSE, hjust = 1, size = 5, color = 'red') +
    geom_text(data = p_values,
              aes(x = max(benchmarking_df$Frequency, na.rm = TRUE) * 1.05,
                  y = Filtering,
                  label = star_adj),
              inherit.aes = FALSE, hjust = 0, size = 5, color = 'darkgreen') +
    theme_bw() +
    ggtitle(paste0(celltype_of_interest)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(size = 8)
    ) +
    scale_fill_manual(values = col_vector[c(19, 30)]) +
    ylab('')
  
  c2 <- ggplot(benchmarking_df, aes(y = Filtering, fill = KRAS)) +
    geom_bar(stat = 'count', position = 'dodge') +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.y = element_blank()
    ) +
    scale_fill_manual(values = col_vector[c(19, 30)]) +
    xlab('Patient') + ylab('')
  
  combined_plot <- ggarrange(
    c1, c2, nrow = 1, align = 'h',
    widths = c(4, 1), common.legend = TRUE, legend = 'bottom'
  )
  
  plot_list[[celltype_of_interest]] <- combined_plot
  
  print(paste0('Plot for ', celltype_of_interest, ' finished!'))
}

final_plot <- ggarrange(
  plotlist = plot_list,
  ncol = 4, nrow = 5, legend = "bottom", common.legend = TRUE
)

final_plot

ggsave2("benchmark_all_celltypes.pdf", plot = final_plot, path = figures_folder,
        width = 40, height = 40, units = "cm", dpi = 600)


#2D scCODA RESULTS---------

sccoda_df <- data.frame(
  Covariate = rep("KRAS_statusT.mut", 17),
  CellType = c(
    "Cycling", "Cytotoxic PD1 high", "Cytotoxic TEMRA", "IFN-induced",
    "PD1 high", "T reg CD25 high", "T reg CD25 low", "Tcm/Naive",
    "Tem Tfh like CXCR5 high CXCL13 high", "Tem Tfh like CXCR5 low CXCL13 high",
    "Tem Th1 like GZMA high", "Tem Th1 like GZMA low", "Tem Th1 like GZMK high",
    "Tem Th17 like CCR6 high IL17 high", "Tem Th17 like CCR6 high IL17 low",
    "Trm Tfh like CXCR5 high CXCL13 low", "Trm Th1 like TNF high"
  ),
  FinalParameter = c(0.000, 0.550, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                     0.527, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000),
  ExpectedSample = c(34.343, 30.784, 24.542, 27.588, 29.003, 131.419, 191.597, 102.145,
                     55.038, 37.652, 74.768, 240.181, 63.713, 21.942, 45.259, 17.226, 57.477),
  log2FC = c(-0.044, 0.749, -0.044, -0.044, -0.044, -0.044, -0.044, -0.044,
             0.716, -0.044, -0.044, -0.044, -0.044, -0.044, -0.044, -0.044, -0.044)
)

sccoda_df$CellType <- factor(sccoda_df$CellType, levels = levels(integrated@meta.data$CelltypeAnnotationLevel2))

ggplot(sccoda_df, aes(x = CellType, y = log2FC, fill = FinalParameter)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "gray80", high = "darkred") +
  labs(x = "", y = "log2 Fold Change", fill = "PIP") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.y = element_blank(),
  ) +
  ylim(-0.05, 0.8)

ggsave2("2D_scCODA_logFC.pdf", path = figures_folder,
        width = 14, height = 10, units = "cm", dpi = 600)

ggplot(sccoda_df, aes(x = CellType, y = log2FC)) +
  geom_segment(aes(xend = CellType, y = 0, yend = log2FC), color = "gray70") +
  geom_point(aes(color = FinalParameter, size = FinalParameter)) +
  scale_color_gradient(low = "gray", high = "red") +
  scale_size(range = c(1, 6)) +
  coord_flip() +
  labs(x = "", y = "log2 Fold Change", color = "PIP", size = "PIP") +
  theme_minimal()

ggplot(sccoda_df, aes(x = log2FC, y = FinalParameter, label = CellType)) +
  geom_point(color = "black") +
  geom_text(data = subset(sccoda_df, FinalParameter > 0.5),
            hjust = -0.1, vjust = 0.5, size = 3) +
  labs(x = "log2 Fold Change", y = "Posterior Inclusion Probability") +
  theme_minimal()

set.seed(123)
sccoda_df_jittered <- sccoda_df[rep(1:nrow(sccoda_df), each = 5), ]
sccoda_df_jittered$log2FC_jitter <- sccoda_df_jittered$log2FC + rnorm(nrow(sccoda_df_jittered), 0, 0.02)

ggplot(sccoda_df_jittered, aes(y = CellType, x = log2FC_jitter)) +
  geom_jitter(height = 0.2, width = 0.05, color = "gray50", size = 1) +
  geom_point(data = sccoda_df, aes(x = log2FC, y = CellType), color = "red", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "log2 Fold Change", y = "") +
  theme_minimal()

# Potential figure composition_kraspre_vs_post ---------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))


figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras"

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_17Apr25.rds')

integrated <- subset(integrated, subset = CelltypeAnnotationLow %in% setdiff(unique(integrated@meta.data$CelltypeAnnotationLow), c('Artefact', 'Unknown')))

integrated@meta.data$Treatment <- ifelse(is.na(integrated@meta.data$Treatment), 'pre', integrated@meta.data$Treatment)


celltype_composition <- integrated@meta.data %>%
  filter(KRAS == 'mut') %>%
  filter(!CelltypeAnnotationLow %in% c('Unknown', 'Artefact')) %>%
  group_by(Treatment, Patient, CelltypeAnnotationLow) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Total = sum(Count)) %>%
  mutate(Frequency = (Count / Total)) %>%
  filter(Frequency < 0.7) %>%
  filter(Total > 300) %>%
  ungroup() %>%
  arrange(Patient, desc(Frequency))

patients_per_celltypes <- celltype_composition %>% group_by(CelltypeAnnotationLow, Treatment) %>% summarise(NumPatients = n_distinct(Patient), .groups = "drop")

y_limit <- ceiling(max(patients_per_celltypes$NumPatients) / 5) * 5

c1 <- ggplot(data = celltype_composition, aes(x = CelltypeAnnotationLow, y = Frequency, fill = Treatment))+
  geom_point(position = position_jitterdodge(), size = 1.5) +
  geom_boxplot(outliers = FALSE) +
  xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = col_vector[c(10,32)])

c2 <- ggplot(data = celltype_composition, aes(x = CelltypeAnnotationLow, fill = Treatment)) +
  geom_bar(stat = 'count', position = 'dodge') +
  xlab('') +
  ylab('Count') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_y_continuous(
    limits = c(0, y_limit),
    breaks = seq(0, y_limit, by = 5)) +
  scale_fill_manual(values = col_vector[c(10,32)])

ggpubr::ggarrange(c2, c1, nrow = 2, align = 'v', common.legend = TRUE, legend = 'right', heights = c(1, 4), vjust = 0.1)

ggsave2("celltype_composition_krasPreVsPost.pdf", path = figures_folder, width = 20, height = 14, units = "cm", dpi = 300)


wilcox_results <- celltype_composition %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  group_by(CelltypeAnnotationLow) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Frequency ~ Treatment)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)


# 2D and pseudoS5C Milo --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(miloR)
library(SingleCellExperiment)

data_folder <- '/projects/sle_jul_23_gabibov/luad/milo/data/'
figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras"


integrated_kras_milo <- readRDS(paste0(data_folder, 'integrated_milo_final.rds'))
da_results <- read.csv(paste0(data_folder, 'da_milo_annotated.csv'))


ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)


da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel1")
da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel2")

da_results$CelltypeAnnotationLevel1 <- factor(da_results$CelltypeAnnotationLevel1, levels = names(subsets_cols_level1))
da_results$CelltypeAnnotationLevel2 <- factor(da_results$CelltypeAnnotationLevel2, levels = setdiff(names(subsets_cols_level2), 'Unknown'))

plotDAbeeswarm(da_results, group.by = "CelltypeAnnotationLevel2")

cowplot::ggsave2('da_plot_final.pdf', path = figures_folder, height = 20, width = 20, units = 'cm', dpi = 300)


#----------------------------------------------------------------FILTER NHOODS BY FREQUENCY IN CLUSTERS

nhood_freq <- 0.95

fdr_threshold <- 0.05

da_results_filtered <- filter(da_results, CelltypeAnnotationLevel2_fraction >= nhood_freq)

da_medians <- da_results_filtered %>%
  group_by(CelltypeAnnotationLevel2) %>%
  summarise(median_logFC = median(logFC, na.rm = TRUE)) %>%
  mutate(group_by = CelltypeAnnotationLevel2)

da_medians$group_by <- factor(da_medians$group_by, levels = unique(da_results_filtered$CelltypeAnnotationLevel2))

da_medians$pos_y <- seq_along(da_medians$group_by)

plotDAbeeswarm(da_results_filtered, group.by = "CelltypeAnnotationLevel2", alpha = fdr_threshold) +
  ggplot2::geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
  geom_point(data = da_medians,
             aes(y = median_logFC, x = pos_y),
             inherit.aes = FALSE,
             color = "red", size = 3) +
  geom_text(data = da_medians,
            aes(y = -4.1, x = pos_y, label = round(median_logFC, 2)),
            inherit.aes = FALSE,
            size = 5,
            hjust = 1) +
  scale_y_continuous(
    limits = c(-5, 4),
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    minor_breaks = NULL
  ) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

cowplot::ggsave2('2D_da_plot_final_filtered.pdf', path = figures_folder, height = 20, width = 25, units = 'cm', dpi = 300)

#----------------------------------------------------------------GOODNESS-OF-FIT TESTS

# The idea here is to test that distribution of number of statistically significant nhoods in KRAS-mut (positive) and KRAS-WT (negative) differs from equal 0.5:0.5

metrics_summary <- da_results_filtered %>%
  mutate(
    is_signif = SpatialFDR < fdr_threshold,
    logFC_dir = case_when(
      is_signif & logFC > 0 ~ "signif_pos",
      is_signif & logFC < 0 ~ "signif_neg",
      TRUE ~ "other"
    )
  ) %>%
  group_by(CelltypeAnnotationLevel2) %>%
  summarise(
    total_nhoods = n(),
    n_significant = sum(is_signif),
    freq_significant = mean(is_signif),
    n_nonsignificant = sum(!is_signif),
    freq_nonsignificant = mean(!is_signif),
    n_significant_logFC_pos = sum(logFC_dir == "signif_pos"),
    n_significant_logFC_neg = sum(logFC_dir == "signif_neg"),
    freq_logFC_pos_total = n_significant_logFC_pos / total_nhoods,
    freq_logFC_neg_total = n_significant_logFC_neg / total_nhoods,
    freq_logFC_pos_within_signif = n_significant_logFC_pos / n_significant,
    freq_logFC_neg_within_signif = n_significant_logFC_neg / n_significant,
    
    # Chi-square goodness-of-fit test
    chisq_gof_p = purrr::map2_dbl(n_significant_logFC_pos, n_significant_logFC_neg, ~ {
      obs <- c(.x, .y)
      if (sum(obs) > 0) {
        chisq.test(obs, p = c(0.5, 0.5))$p.value
      } else {
        NA_real_
      }
    }),
    
    # Fisher's exact test version of the same logic (observed vs expected counts)
    fisher_gof_p = purrr::map2_dbl(n_significant_logFC_pos, n_significant_logFC_neg, ~ {
      total <- .x + .y
      if (total > 0) {
        expected <- total / 2
        tbl <- matrix(c(.x, .y, expected, expected), nrow = 2)
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    })
  ) %>%
  # Apply multiple testing correction
  mutate(
    chisq_gof_p_adj_bh = p.adjust(chisq_gof_p, method = "BH"),
    fisher_gof_p_adj_bh = p.adjust(fisher_gof_p, method = "BH"),
    chisq_gof_p_adj_bonf = p.adjust(chisq_gof_p, method = "bonferroni"),
    fisher_gof_p_adj_bonf = p.adjust(fisher_gof_p, method = "bonferroni")
  )

#----------------------------------------------------------------PLOT NHOODS ON KNN GRAPH

integrated_kras_milo <- buildNhoodGraph(integrated_kras_milo)

combined_plot <- DimPlot(integrated_kras, group.by = 'CelltypeAnnotationLevel2', raster = FALSE, cols = subsets_cols_level2, pt.size = 2) + 
  plotNhoodGraphDA(integrated_kras_milo, da_results, alpha=0.05, size_range=c(2,8), highlight.da = 0) +
  plot_layout(guides = "collect") #+
#  geom_point_rast(size = 0.3, raster.dpi = 300)

cowplot::ggsave2(filename = file.path(figures_folder, 'umaps_final.pdf'),
                 plot = combined_plot,
                 height = 25, width = 50, units = 'cm', dpi = 300)

sig_nhoods <- as.character(da_results[da_results$SpatialFDR < fdr_threshold, "Nhood"])

combined_plot2 <- DimPlot(integrated_kras, group.by = 'CelltypeAnnotationLevel2', raster = FALSE, cols = subsets_cols_level2) + 
  plotNhoodGraphDA(integrated_kras_milo, da_results, subset.nhoods = sig_nhoods, alpha=0.05) +
  plot_layout(guides = "collect")

cowplot::ggsave2(filename = file.path(figures_folder, 'pseudoS5C_umaps_final_filtrered.pdf'),
                 plot = combined_plot2,
                 height = 30, width = 60, units = 'cm', dpi = 300)


#S6A S6B S7 Tregs -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(kableExtra)
library(cowplot)

source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'


figures_folder <- '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras'
data_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/data/'
rds_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/rds/'

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')


integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel2

#-----------------------------------------------------------------SUBSET TREGS

tregs <- subset(integrated, subset = CelltypeAnnotationLow %in% c('T reg CD25 high' , 'T reg CD25 low', 'IFN-induced'))

tregs[["RNA"]] <- split(tregs[["RNA"]], f = tregs$orig.ident)

preprocessing <- function(obj){
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- ScaleData(obj)
  variable_genes <- VariableFeatures(object = obj)
  #genes_to_remove <- "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*|^IGK*|^IGL*|^IGH*|TRBC1|TRBC2|TRDC|TRGC1|TRGC2|^MT"
  genes_to_remove <- "^TRBV|^TRBD|^TRBJ|^TRDV|^TRDD|^TRDJ|^TRAV|^TRAJ|^TRGV|^TRGJ|^IGK|^IGL|^IGH|TRBC1|TRBC2|TRDC|TRGC1|TRGC2|^MT"
  variable_genes <- variable_genes[-grep(genes_to_remove, variable_genes, perl = TRUE)]
  variable_genes <- head(variable_genes, 2000)
  print(variable_genes)
  obj <- RunPCA(obj, features = variable_genes)
  return (obj)
}

tregs <- preprocessing(tregs)


integrated_tregs <- IntegrateLayers(
  object = tregs, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca"
)

integrated_tregs[["RNA"]] <- JoinLayers(integrated_tregs[["RNA"]])

ElbowPlot(integrated_tregs, ndims = 30, reduction = "pca")

integrated_tregs <- FindNeighbors(integrated_tregs, reduction = "integrated.rpca", dims = 1:30)
integrated_tregs <- FindClusters(integrated_tregs, resolution = 0.5)

integrated_tregs <- RunUMAP(integrated_tregs, reduction = "integrated.rpca", dims = 1:30)
# Skip here
integrated_tregs <- readRDS(paste0(rds_folder, 'integrated_tregs.rds'))
integrated_tregs <- subset(integrated_tregs, idents = '10', invert = T) # removing IFN induced

DimPlot(integrated_tregs, group.by='seurat_clusters', cols=col_vector, raster = FALSE, label = T) +
  xlab('UMAP1') +
  ylab('UMAP2')
ggsave2("S6A_clusters_umap.pdf", path = figures_folder, width = 15, height = 10, units = "cm", dpi = 600)

#-----------------------------------------------------------------SAVE RDS

#saveRDS(integrated_tregs, paste0(rds_folder, 'integrated_tregs.rds')) ## For some reason (probably due to updates in Seurat or in the original object), clustering work differently since original generation. Need to look out for this

#-----------------------------------------------------------------EXPRESSION OF MARKER GENES
integrated_tregs <- readRDS(paste0(rds_folder, 'integrated_tregs.rds'))
integrated_tregs <- subset(integrated_tregs, idents = '10', invert = T) # removing IFN induced

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

genes <- unique(unlist(gene_signatures))
cell_types <- rep(names(gene_signatures), sapply(gene_signatures, length))
gene_df <- data.frame(Cell_Type = cell_types, Gene = unlist(gene_signatures))

DotPlot(integrated_tregs, features = unique(gene_df$Gene), group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave2("marker_genes.png", path = figures_folder, width = 45, height = 30, units = "cm")
ggsave2("marker_genes.pdf", path = figures_folder, width = 45, height = 30, units = "cm")

#-----------------------------------------------------------------UMAP OF MARKER GENES

umap_features <- c('FOXP3', 'IKZF2', 'IL2RA', 'CCL5', 'GZMK', 'IFNG', 'PDCD1', 'LAG3', 'CTLA4', 'TIGIT')
integrated_tregs <- subset(integrated_tregs, idents = '10', invert = T) # removing IFN induced

p_cluster <- DimPlot(integrated_tregs, group.by = "seurat_clusters", raster = FALSE, cols = col_vector, label = TRUE)

feature_plots <- lapply(umap_features, function(feature) {
  FeaturePlot(integrated_tregs, features = feature, raster = T, order = FALSE)
})

all_plots <- wrap_plots(c(list(p_cluster), feature_plots), ncol = 4)

ggsave2(filename = "treg_combined_umap_plots.pdf",
        path = figures_folder,
        plot = all_plots,
        width = 30,
        height = 20,
        dpi = 300)


integrated_tregs@meta.data$orig.ident_old <- integrated_tregs@meta.data$orig.ident

integrated_tregs@meta.data <- integrated_tregs@meta.data %>%
  mutate(orig.ident_custom = case_when(
    stringr::str_starts(orig.ident, 'P') ~ 'WetLab data',
    TRUE ~ orig.ident
  ))

DimPlot(integrated_tregs, group.by='orig.ident', cols=col_vector, raster = FALSE) +
  facet_wrap(facet = 'orig.ident')
ggsave2("datasets_split_umap_tregs.pdf", path = figures_folder, width = 30, height = 30, units = "cm")

# Plot with subsampling
subsampling_number <- 1500

cells_to_keep <- integrated_tregs@meta.data %>%
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

integrated_tregs_subsampled <- subset(integrated_tregs, cells = cells_to_keep)

DimPlot(integrated_tregs_subsampled, group.by='seurat_clusters', split.by = 'orig.ident_custom', cols=col_vector, raster=FALSE, pt.size = 1.5) +
  facet_wrap(facet = 'orig.ident_custom', ncol = 4) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position = 'bottom')

ggsave2("S7_datasets_split_umap_subsampled1_5k_tregs.pdf", path = figures_folder, width = 45, height = 30, units = "cm", dpi = 600)
# Plot with subsampling

plot_data_unnorm <- plot_batch_fr(integrated_tregs, cluster_column = 'seurat_clusters', batch_column = 'orig.ident', palette = col_vector, normalize = FALSE)
ggsave2("S6B_batch_disribution_over_clusters_treg.pdf", path = figures_folder, width = 18, height = 10, units = "cm")

plot_data_norm <- plot_batch_fr(integrated_tregs, cluster_column = 'seurat_clusters', batch_column = 'orig.ident_custom', palette = col_vector, normalize = TRUE)
ggsave2("S6B_normalized_batch_disribution_over_clusters_treg.pdf", path = figures_folder, width = 18, height = 10, units = "cm")

#-----------------------------------------------------------------DEG HEATMAP TREGS

Idents(integrated_tregs) <- integrated_tregs@meta.data$seurat_clusters

markers <- FindAllMarkers(integrated_tregs,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.58)

#write.csv(markers, paste0(data_folder, 'degs_tregs.csv'))

# markers <- read.csv(paste0(data_folder, 'degs_per_celltype.csv'))

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

subsampling_number <- 500

cells_to_keep500 <- integrated_tregs@meta.data %>%
  tibble::rownames_to_column(var = 'barcodes') %>%
  group_by(seurat_clusters) %>%
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

integrated_tregs_subsampled500 <- subset(integrated_tregs, cells = cells_to_keep500)

Idents(integrated_tregs_subsampled500) <- integrated_tregs_subsampled500@meta.data$seurat_clusters

DoHeatmap(integrated_tregs_subsampled500, features = top10$gene, size = 4, group.colors = col_vector) +
  scale_fill_gradientn(colors = c("purple", "black", "yellow")) +
  theme(axis.text.x = element_blank()) + 
  guides(color = "none", fill = guide_colorbar(title = "Expression"))

cowplot::ggsave2('heatmap5top_tregs.pdf', path = figures_folder, width = 30, height = 40, units = 'cm', dpi = 300)

#3AB  S6C S6D 3C STUDY Type1 Treg CLUSTER------------------------------------------------------

# COMPARE MARKERS OF T REGS AND TH1 ACROSS TH1 AND NEW SUBSET

library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(kableExtra)
library(cowplot)

source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'


figures_folder <- '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras'
data_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/data/'
rds_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/rds/'

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')


integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel2

integrated_tregs <- readRDS(paste0(rds_folder, 'integrated_tregs.rds'))
integrated_tregs <- subset(integrated_tregs, idents = '10', invert = T) # removing IFN induced
cells_type1_treg <- WhichCells(integrated_tregs, idents = '8')
comparison <- subset(integrated, subset = CelltypeAnnotationLevel1 %in% c('Tem Th1 like', 'Trm Th1 like', 'T reg'))
Idents(comparison) <- 'CelltypeAnnotationLevel1'
Idents(comparison, cells = cells_type1_treg) <- 'Type-1 Treg'
comparison$Subset <- Idents(comparison) # to keep needed idents

# 0. Setup: fix cluster order so Type-1 Treg sits next to the Th1s
ident_col <- "Subset"   # <-- your metadata column

order_lvls <- c("Trm Th1 like",
                "Tem Th1 like",
                "Type-1 Treg",
                "T reg")

seu <- comparison   # <-- your object
seu@meta.data[[ident_col]] <- factor(seu@meta.data[[ident_col]],
                                     levels = order_lvls)
Idents(seu) <- ident_col
DefaultAssay(seu) <- "RNA"

# Use log-normalised data for dot plots and module scores.
# If you only ran SCTransform, normalise RNA too:
if (!"data" %in% slotNames(seu[["RNA"]]) ||
    all(dim(GetAssayData(seu, assay = "RNA", slot = "data")) == 0)) {
  seu <- NormalizeData(seu, assay = "RNA")
}

# helper: keep only genes actually present, warn about the rest
present <- function(g) {
  keep <- g[g %in% rownames(seu)]
  miss <- setdiff(g, keep)
  if (length(miss)) message("Not in object, dropped: ", paste(miss, collapse = ", "))
  keep
}



# ---- define genes with their group labels 
gene_groups <- list(
  panel_treg = c(
    "FOXP3" = "Core TF",  "IKZF2" = "Core TF",
    "IKZF4" = "Core TF", "STAT5B" = "Core TF",
    "IL2RA" = "Surface suppressive function", "CTLA4" = "Surface suppressive function", 
    "TIGIT" = "Surface suppressive function",
    "CCR8"  = "Surface suppressive function", "IL1R2" = "Surface suppressive function"
  ),
  panel_eff = c(
    "IFNG" = "Cytokines", "TNF" = "Cytokines",
    "CCL5" = "Chemokines", "CCL4" = "Chemokines",
    "GZMK" = "Cytotoxic granules / licensing", "GZMA" = "Cytotoxic granules / licensing", "GZMH" = "Cytotoxic granules / licensing",
    "NKG7" = "Cytotoxic granules / licensing", "CTSW" = "Cytotoxic granules / licensing",
    "SYTL2" = "Cytotoxic granules / licensing" 
  ),
  panel_tf = c(
    "TBX21" = "Type-1 drivers", "STAT4" = "Type-1 drivers",
    "BHLHE40" = "Type-1 drivers", "EOMES" = "Type-1 drivers", 
    "RUNX3" = "Terminal effector / residency", "ID2" = "Terminal effector / residency",
    "ZEB2"  = "Terminal effector / residency", "KLRG1" = "Terminal effector / residency",
    "HOPX"  = "Terminal effector / residency", "PRDM1" = "Terminal effector / residency",
    "TCF7" = "Stem-like", "LEF1" = "Stem-like",
    "GATA3" = "Other lineages", "RORC" = "Other lineages", "BCL6" = "Other lineages"
  ),
  panel_disc = c(
    "ANXA1" = "Activation", "CD40LG" = "Activation", "NR4A1" = "Activation", "NR4A2" = "Activation",
      "PDCD1"  = "Cis-inhibition",
      "HAVCR2" = "Cis-inhibition",
      "KLRG1"  = "Cis-inhibition", 
      "LAG3"   = "Trans-inhibition",
      "ENTPD1" = "Trans-inhibition",
    "CD69" = "Trafficking", "CCR7" = "Trafficking",
    "SELL" = "Trafficking", "CXCR3" = "Trafficking", 
    "CXCR6" = "Trafficking","CCR5" = "Trafficking", 
    "S1PR1" = "Trafficking"
    )
  )


# drop genes absent from the object, keep order
clean_group <- function(g) g[names(g) %in% rownames(seu)]
gene_groups <- lapply(gene_groups, clean_group)

# ---- faceted dot plot builder 
make_dot_faceted <- function(gmap, title, scale = TRUE,
                             pct_lim = c(0, 100),
                             pct_breaks = c(0, 25, 50, 75, 100)) {
  genes <- names(gmap)
  p <- DotPlot(seu, features = genes, scale = scale,
               dot.scale = 6, cols = c("lightgrey", "darkred"))
  p$data$gene_group <- factor(gmap[as.character(p$data$features.plot)],
                              levels = unique(gmap))
  
  p +
    scale_size(
      limits = pct_lim, range = c(0, 6), breaks = pct_breaks,
      name = "Percent\nexpressed"
    )  +
    facet_grid(~ gene_group, scales = "free_x", space = "free_x") +
    labs(title = title) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.title   = element_blank(),
      panel.grid   = element_line(linewidth = 0.2),
      panel.spacing.x = unit(2, "pt"),
      strip.background = element_rect(fill = "grey92", colour = "grey60"),
      strip.text   = element_text(size = 7.5, lineheight = 0.9,
                                  margin = margin(2, 2, 2, 2)),
      plot.title   = element_text(size = 10, face = "bold"),
      legend.key.size = unit(0.35, "cm")
    )
}


d1 <- make_dot_faceted(gene_groups$panel_treg, "T Regulatory identity")
d2 <- make_dot_faceted(gene_groups$panel_eff,  "Type-1 / cytotoxic effector program")
d3 <- make_dot_faceted(gene_groups$panel_tf,   "Transcription factors")
d4 <- make_dot_faceted(gene_groups$panel_disc, "Activation / Cis-trans inhibition/ trafficking")

dot_fig <- (d1 / d2 / d3 / d4) +
  plot_layout(guides = "collect", heights = c(1, 1, 1.05, 1)) &
  theme(legend.position = "right")

ggsave(paste0(figures_folder, "/3A_fig_dotplot_grouped.pdf"), dot_fig, width = 9, height = 8,
       useDingbats = FALSE)

 # Gene modules
library(grid)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rlang)
library(rstatix)
library(ggpubr)

# 1. Module gene sets
mods <- list(
  Treg_score  = present(c("FOXP3","IL2RA","IKZF2", "CCR8",
                          "IL1R2")),
  Cytotox_score = present(c("GZMA","GZMK","GZMH","GZMB","PRF1",
                            "NKG7","CTSW","FASLG","GNLY","KLRD1")),
  Type1_score = present(c("IFNG","TBX21","EOMES","STAT4",
                          "CCL4","CCL5","IFNGR1","BHLHE40"))
)

set.seed(1)
seu <- AddModuleScore(seu, features = mods, name = "MS_",
                      assay = "RNA", ctrl = 50)

ms_cols <- paste0("MS_", seq_along(mods))
colnames(seu@meta.data)[match(ms_cols, colnames(seu@meta.data))] <- names(mods)
# 2. Display order on the violin x axis (reversed vs the dot plots)
vln_levels <- rev(order_lvls)

md <- seu@meta.data |>
  dplyr::select(group = !!ident_col, dplyr::all_of(names(mods))) |>
  dplyr::mutate(group = factor(as.character(group), levels = vln_levels))

# 3. Per-donor pseudobulk: donors are the unit of analysis

donor_col <- "Patient"   # <-- your donor / sample column

pb <- seu@meta.data |>
  dplyr::group_by(donor = .data[[donor_col]],
                  group = .data[[ident_col]]) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   dplyr::across(dplyr::all_of(names(mods)), mean),
                   .groups = "drop") |>
  dplyr::filter(n_cells >= 10) |>            # drop donor-subset pairs with too few cells
  dplyr::mutate(group = factor(as.character(group), levels = vln_levels))

# donors contributing all four subsets -> paired test is available
complete_donors <- pb |>
  dplyr::count(donor) |>
  dplyr::filter(n == length(order_lvls)) |>
  dplyr::pull(donor)

use_paired <- length(complete_donors) >= 3
pb_test <- if (use_paired) dplyr::filter(pb, donor %in% complete_donors) else pb

message("Donors passing filter: ", dplyr::n_distinct(pb$donor),
        " | complete across all subsets: ", length(complete_donors),
        " | paired test: ", use_paired)

# 4. Per-donor Wilcoxon vs Type-1 Treg, one BH correction over all scores
stat_pb <- lapply(names(mods), function(sc) {
  pb_test |>
    dplyr::select(donor, group, val = !!sym(sc)) |>
    dplyr::arrange(group, donor) |>
    rstatix::wilcox_test(val ~ group, ref.group = "Type-1 Treg",
                         paired = use_paired) |>
    dplyr::mutate(score = sc)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p.adj  = p.adjust(p, method = "BH"),
    signif = dplyr::case_when(
      p.adj <= 1e-4 ~ "****",
      p.adj <= 1e-3 ~ "***",
      p.adj <= 1e-2 ~ "**",
      p.adj <= 0.05 ~ "*",
      TRUE          ~ "ns"
    )
  )

print(stat_pb)
write.csv(stat_pb, paste0(figures_folder,"/module_score_stats_perdonor.csv"), row.names = FALSE)

# 5. Violin panel: cell-level violins, per-donor p-values as brackets
make_vln <- function(sc, title, hide_ns = FALSE) {
  df <- md |> dplyr::select(group, val = !!sym(sc))
  
  st <- stat_pb |> dplyr::filter(score == sc)
  if (hide_ns) st <- dplyr::filter(st, signif != "ns")
  
  yr <- diff(range(df$val, na.rm = TRUE))
  st$xmin       <- match(st$group1, vln_levels)
  st$xmax       <- match(st$group2, vln_levels)
  st$y.position <- max(df$val, na.rm = TRUE) + yr * 0.09 * seq_len(nrow(st))
  st$xmid       <- (st$xmin + st$xmax) / 2
  
  p <- VlnPlot(seu, features = sc, pt.size = 0, sort = FALSE) + NoLegend()
  
  p <- p +
    stat_summary(data = p$data,
                 aes(x = ident, y = .data[[sc]]),
                 inherit.aes = FALSE,
                 fun = median, geom = "crossbar",
                 width = 0.5, linewidth = 0.3, colour = "black") +
    scale_x_discrete(limits = vln_levels) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = title, y = "Module score") 
  
  if (nrow(st)) {
    p <- p +
      geom_segment(data = st,
                   aes(x = xmin, xend = xmax,
                       y = y.position, yend = y.position),
                   inherit.aes = FALSE,
                   linewidth = 0.3, colour = "black") +
      geom_text(data = st,
                aes(x = xmid, y = y.position, label = signif),
                inherit.aes = FALSE,
                size = 3.2, vjust = -0.25, colour = "black")
  }
  p
}

v1 <- make_vln("Treg_score",    "Treg module")
v2 <- make_vln("Cytotox_score", "Cytotoxicity module")
v3 <- make_vln("Type1_score",   "Type-1 module")
#no_x <- theme(axis.text.x = element_blank(),
#              axis.ticks.x = element_blank())

vln_fig <- (v1 + no_x) / (v2 + no_x) / v3 + plot_layout(nrow = 1) &
  theme(legend.position = "none")
ggsave(paste0(figures_folder, "/3B_vln_modules_type1_treg.pdf"), vln_fig, width = 9, height = 3.2,
       useDingbats = FALSE)

# Differential expression on Tregs - IFN induced are removed previosuly
th1_treg_markers <- FindMarkers(integrated_tregs, ident.1 = '8', only.pos = F, logFC.threshold = 0.2,
                                )
th1_treg_markers <- filter(th1_treg_markers,  p_val_adj < 0.01, abs(avg_log2FC) > 0.2)
th1_treg_markers$gene <- rownames(th1_treg_markers)
th1_treg_markers <- th1_treg_markers[, c("gene", setdiff(names(th1_treg_markers), "gene"))]  # gene first
writexl::write_xlsx(th1_treg_markers, paste0(figures_folder, "/type1_treg_markers.xlsx"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                           rownames(qual_col_pals)))

col_vector[4] <- '#FDDA0D'
integrated_tregs <- RenameIdents(integrated_tregs, '8' = 'Type-1 Treg')


umap_features <- c('FOXP3', 'IKZF2', 'IL2RA', 'CCL5', 'GZMK', 'IFNG', 'EOMES', 'ANXA1', 'CD40LG', 'HOPX', 'ITGA1')

p_cluster <- DimPlot(integrated_tregs, raster = F, cols = col_vector, label = T)

feature_plots <- lapply(umap_features, function(feature) {
  FeaturePlot(integrated_tregs, features = feature, raster = F, order = T)
})

all_plots <- wrap_plots(c(list(p_cluster), feature_plots), ncol = 3)

ggsave2(filename = "treg_combined_umap_plots.png",
        path = figures_folder,
        plot = all_plots,
        width = 15,
        height = 15,
        dpi = 300)

# S6? Treg Plasticity of type1 treg --------------------------------------------
# Plasticity for Treg subset intra-Treg
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(forcats)
rds_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/rds/'
figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed"

intgr <- readRDS(paste0(rds_folder, 'integrated_tregs.rds'))
intgr <- subset(intgr, idents = '10', invert = T) # removing IFN-induced
intgr@meta.data$patient_id <- intgr@meta.data$Patient
intgr@meta.data$Subset <- intgr@meta.data$seurat_clusters
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$TCR2nt[intgr$TCR2nt == "NA"] <- NA
tcr_obj_cells <- dplyr::select(intgr[[]], TCR2nt) %>% drop_na()
intgr <- subset(intgr, cells = row.names(tcr_obj_cells))

intgr$donor_TCR2nt <- paste0(intgr$patient_id, "_", intgr$TCR2nt) # Count only intra-patient plasticity
# No less than N clonotypes in each Subset in each Patient
N = 50
intgr$Barcode <- row.names(intgr[[]])
intgr_filtered_table <- intgr[[]] %>%
  group_by(Subset, patient_id) %>%
  filter(n_distinct(TCR2nt) > N) %>%
  ungroup()
intgr <- subset(intgr, cells = intgr_filtered_table$Barcode)
K <- 2

clones_per_observation <- intgr@meta.data %>% group_by(Subset, patient_id) %>% summarize(n = n())
patients_number <- intgr@meta.data %>% group_by(Subset) %>% summarize(n = n_distinct(patient_id))

relevant_celltypes <- patients_number %>% filter(n >= K) %>% pull(Subset)

intgr <- subset(intgr, subset = seurat_clusters %in% relevant_celltypes)


#write.csv(patients_number, paste0(data_folder, 'plasticityWT_Combined_patientPerCelltype_', N, 'cells.csv'))

# Plotting D metric

plotPlasticity <- function(obj, saveMat = TRUE, pathMat = NULL){
  cell_type.tcr.list <- list()
  for (cell_type in unique(obj@meta.data$Subset)) {
    i <- obj@meta.data[obj@meta.data$Subset == cell_type, ]$donor_TCR2nt
    i <- i[!is.na(i)]
    cell_type.tcr.list[[cell_type]] <- i
  }
  D.old <- function(i, j, C=1)
    #i - вектор последовательностей TRB клеток кластера i
    #j - вектор последовательностей TRB клеток кластера j
    #C - константа (по умолчанию - 1)
  {
    result <- length(intersect(i, j)) * C/(as.numeric(n_distinct(i)) * as.numeric(n_distinct(j)))
    return(result)
  }
  graph_columns <- c("clusterA", "clusterB", "D", "clusterA_size")
  graph_table <- data.frame(matrix(nrow = 0, ncol = length(graph_columns)))
  colnames(graph_table) <- graph_columns
  for (i in names(cell_type.tcr.list)) {
    s <- length(cell_type.tcr.list[[i]])
    for (j in names(cell_type.tcr.list)) {
      d <- round(D.old(cell_type.tcr.list[[i]], cell_type.tcr.list[[j]]), 10)
      graph_table[nrow(graph_table)+1, ] <- c(i, j, d, s)
    }
  }
  graph_table$D <- as.numeric(graph_table$D)
  graph_table$logD <- log2(1+graph_table$D)
  # logD metric
  mat <- xtabs(logD~clusterA+clusterB, data = graph_table)
  clust_order <- c('0',
                   '1',
                   '2',
                   '3',
                   '4',
                   '5',
                   '6',
                   '7',
                   '8',
                   '9',
                   '10',
                   '11',
                   '12')
  
  cluster_cols <- c(
    '0' = '#1B9E77',                       # Teal green
    '1' = '#66C2A5',                      # Lighter teal
    '2' = '#01665E',                        # Dark cyan
    '3' = '#1F78B4',          # Deep blue
    '4' = '#6A3D9A',    # Dark purple
    '5' = '#B15928',          # Brownish-orange
    '6' = '#FF7F00',               # Bright orange
    '7' = '#FB9A99',                      # Light red/pink
    '8' = '#8DD3C7',    # Aqua green
    '9' = '#FDDA0D',         # Yellow
    '10' = '#984EA3',             # Dark purple
    '11' = '#999999',                # Neutral gray
    '12' = '#7D6231'                 # Dark brown
  )
  
  subsets <- colnames(mat)[match(clust_order, colnames(mat), nomatch = 0L)]
  mat <- mat[subsets, subsets]
  subsets <- colnames(mat)
  diag(mat) <- NA
  
  if (saveMat){
    saveRDS(mat, file = pathMat)
  }
  
  col_fun = colorRamp2(c(0, quantile(mat, 0.99, na.rm = T)), c("white", "red"))
  p1 <- Heatmap(mat,
                column_title = expression(" "),
                column_order = subsets,
                column_names_side = "top",
                column_names_rot = 45,
                row_order = subsets,
                show_row_names = T,
                row_names_side = "left",
                width = unit(16.6, "cm"),
                height = unit(16.6, "cm"),
                rect_gp = gpar(col = "gray", lwd = 1),
                col = col_fun,
                heatmap_legend_param = list(title = expression("log D-metric")),
                top_annotation = HeatmapAnnotation(clusters = subsets,
                                                   col = list(clusters = cluster_cols),
                                                   show_annotation_name = F,
                                                   show_legend = F),
                left_annotation = rowAnnotation(clusters = subsets,
                                                col = list(clusters = cluster_cols),
                                                show_legend = F,
                                                show_annotation_name = F))
  plot(p1)
}

plot_D_metric <- plotPlasticity(intgr, saveMat = F, pathMat = '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/')
plot(plot_D_metric)

pdf(file=paste0(figures_folder, paste0("/intra_treg_ifn_induced", N, "clonesF_", K, "_plasticity.pdf")), width=11,height=11)
draw(plot_D_metric)
dev.off()


# Plascitity of type1 Treg with the entire dataset 

library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(forcats)
source('/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/luad_figures_tables/colors_fixed.R')
figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed"

rds_folder <- '/projects/sle_jul_23_gabibov/luad/subclustering/rds/'
treg <- readRDS(paste0(rds_folder, 'integrated_tregs.rds'))
type1_treg <- WhichCells(subset(treg, idents = '8'))
intgr <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')
#intgr$CelltypeAnnotationLow <- intgr$CelltypeAnnotationLevel2
intgr <- subset(intgr, subset = CelltypeAnnotationLevel2 %in% setdiff(unique(intgr@meta.data$CelltypeAnnotationLevel2), c('Unknown')))
Idents(intgr) <- intgr$CelltypeAnnotationLevel2
Idents(intgr, cells = type1_treg) <- 'Type 1-like Treg'
intgr$Subset <- Idents(intgr)

intgr@meta.data$patient_id <- intgr@meta.data$Patient
#intgr@meta.data$Subset <- intgr@meta.data$CelltypeAnnotationLevel2
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$TCR2nt[intgr$TCR2nt == "NA"] <- NA
tcr_obj_cells <- dplyr::select(intgr[[]], TCR2nt) %>% drop_na()
intgr <- subset(intgr, cells = row.names(tcr_obj_cells))

intgr$donor_TCR2nt <- paste0(intgr$patient_id, "_", intgr$TCR2nt) # Count only intra-patient plasticity
# No less than N clonotypes in each Subset in each Patient
N = 50
intgr$Barcode <- row.names(intgr[[]])
intgr_filtered_table <- intgr[[]] %>%
  group_by(Subset, patient_id) %>%
  filter(n_distinct(TCR2nt) > N) %>%
  ungroup()
intgr <- subset(intgr, cells = intgr_filtered_table$Barcode)
K <- 2

clones_per_observation <- intgr@meta.data %>% group_by(Subset, patient_id) %>% summarize(n = n())
patients_number <- intgr@meta.data %>% group_by(Subset) %>% summarize(n = n_distinct(patient_id))

relevant_celltypes <- patients_number %>% filter(n >= K) %>% pull(Subset)

intgr_filtered <- subset(intgr, subset = Subset %in% relevant_celltypes)
Idents(intgr_filtered) <- 'Subset'
intgr_filtered$Subset <- Idents(intgr_filtered)
subsets_cols_level2 <- c(
  'Tcm/Naive'                                     = '#f7b6d2',  # Sky blue
  
  # Th17-like → Greens
  'Tem Th17 like CCR6 high IL17 low'                   = '#6BBE44',  # Light green
  'Tem Th17 like CCR6 high IL17 high'                  = '#228B22',  # Forest green
  
  # Th1-like → Oranges/Reds
  'Trm Th1 like TNF high'                         = '#8DD3C7',
  'Tem Th1 like GZMA low'                         = '#F07F44',
  'Tem Th1 like GZMA high'                        = '#d6616b',
  'Tem Th1 like GZMK high'                        = '#C94D00',
  
  # Cytotoxic
  #  'Cytotoxic TEMRA'                               = '#E603DC',
  'Cytotoxic PD1 high'                            = '#E31A1C',
  
  # Tfh-like → Blues/Teals
  'Trm Tfh like CXCR5 high CXCL13 low'            = '#b8a35a',
  'Tem Tfh like CXCR5 low CXCL13 high'            = '#33AACC',
  'Tem Tfh like CXCR5 high CXCL13 high'           = '#005F87',
  
  # PD1 high
  'PD1 high'                                      = '#FDDA0D',
  
  # T regs
  'T reg CD25 low'                                = '#984EA3',
  'T reg CD25 high'                               = '#BEAED4',
  'Type 1-like Treg'                              = 'darkred',
  
  # IFN-induced
  'IFN-induced'                                   = '#999999',
  
  # Cycling
  'Cycling'                                       = '#E6AB02'
  
  # Unknown
  #  'Unknown'                                       = '#CCCCCC'
)
#write.csv(patients_number, paste0(data_folder, 'plasticityWT_Combined_patientPerCelltype_', N, 'cells.csv'))

# Plotting D metric

plotPlasticity <- function(obj, saveMat = TRUE, pathMat = NULL){
  cell_type.tcr.list <- list()
  for (cell_type in unique(obj@meta.data$Subset)) {
    i <- obj@meta.data[obj@meta.data$Subset == cell_type, ]$donor_TCR2nt
    i <- i[!is.na(i)]
    cell_type.tcr.list[[cell_type]] <- i
  }
  D.old <- function(i, j, C=1)
    #i - вектор последовательностей TRB клеток кластера i
    #j - вектор последовательностей TRB клеток кластера j
    #C - константа (по умолчанию - 1)
  {
    result <- length(intersect(i, j)) * C/(as.numeric(n_distinct(i)) * as.numeric(n_distinct(j)))
    return(result)
  }
  graph_columns <- c("clusterA", "clusterB", "D", "clusterA_size")
  graph_table <- data.frame(matrix(nrow = 0, ncol = length(graph_columns)))
  colnames(graph_table) <- graph_columns
  for (i in names(cell_type.tcr.list)) {
    s <- length(cell_type.tcr.list[[i]])
    for (j in names(cell_type.tcr.list)) {
      d <- round(D.old(cell_type.tcr.list[[i]], cell_type.tcr.list[[j]]), 10)
      graph_table[nrow(graph_table)+1, ] <- c(i, j, d, s)
    }
  }
  graph_table$D <- as.numeric(graph_table$D)
  graph_table$logD <- log2(1+graph_table$D)
  # logD metric
  mat <- xtabs(logD~clusterA+clusterB, data = graph_table)
  clust_order <- names(subsets_cols_level2)
  
  cluster_cols <- subsets_cols_level2                 
  
  
  subsets <- colnames(mat)[match(clust_order, colnames(mat), nomatch = 0L)]
  mat <- mat[subsets, subsets]
  subsets <- colnames(mat)
  diag(mat) <- NA
  
  if (saveMat){
    saveRDS(mat, file = pathMat)
  }
  
  col_fun = colorRamp2(c(0, quantile(mat, 0.99, na.rm = T)), c("white", "red"))
  p1 <- Heatmap(mat,
                column_title = expression(" "),
                column_order = subsets,
                column_names_side = "top",
                column_names_rot = 45,
                row_order = subsets,
                show_row_names = T,
                row_names_side = "left",
                width = unit(16.6, "cm"),
                height = unit(16.6, "cm"),
                rect_gp = gpar(col = "gray", lwd = 1),
                col = col_fun,
                heatmap_legend_param = list(title = expression("log D-metric")),
                top_annotation = HeatmapAnnotation(clusters = subsets,
                                                   col = list(clusters = cluster_cols),
                                                   show_annotation_name = F,
                                                   show_legend = F),
                left_annotation = rowAnnotation(clusters = subsets,
                                                col = list(clusters = cluster_cols),
                                                show_legend = F,
                                                show_annotation_name = F))
  plot(p1)
}

plot_D_metric <- plotPlasticity(intgr, saveMat = F, pathMat = '/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/')
plot(plot_D_metric)

pdf(file=paste0(figures_folder, paste0("/whole_dataset_treg_type1", N, "clonesF_", K, "plasticity_level2.pdf")), width=11,height=11)
draw(plot_D_metric)
dev.off()


 # 3C KRAS mut vs WT frequency 
cluster_of_interest <- 8

batch_patient_stats <- integrated_tregs@meta.data %>% 
  group_by(seurat_clusters, orig.ident_custom, Patient, KRAS) %>%
  summarize(n = n())

batch_patient_cluster_gzma <- batch_patient_stats %>% filter(seurat_clusters == cluster_of_interest) %>% filter(n >= 10)

patients_cluster_gzma <- batch_patient_cluster_gzma %>% pull(Patient)

patients_cluster_gzma_total_cells <- integrated@meta.data %>% filter(Patient %in% patients_cluster_gzma) %>% group_by(Patient) %>% summarize(total = n())

patients_cluster_gzma_total_tregs <- integrated_tregs@meta.data %>% filter(Patient %in% patients_cluster_gzma) %>% group_by(Patient) %>% summarize(total_treg = n())

batch_patient_cluster_gzma <- batch_patient_cluster_gzma %>%
  left_join(patients_cluster_gzma_total_cells, by = 'Patient') %>%
  left_join(patients_cluster_gzma_total_tregs, by = 'Patient')

df_summary <- batch_patient_cluster_gzma %>%
  group_by(orig.ident_custom) %>%
  summarise(
    unique_patients = n_distinct(Patient),
    mean_cells_per_patient = mean(n, na.rm = TRUE)
  ) %>%
  arrange(desc(unique_patients))

ggplot(df_summary, aes(y = reorder(orig.ident_custom, unique_patients), x = unique_patients, fill = mean_cells_per_patient)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradientn(colors = c("blue", "#FDDA0D", "red"), name = "Mean cells\nper patient") +
  labs(x = "Number of patients", y = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 'dashed'))

cowplot::ggsave2('S6C_patients_treg_type1.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 600)

# Check if the amount of this specific cells depend on the amount of total cells in Patient. In my point of view, this means that is can't find this cells 
# if the amount of total cells in Patient is low if there is an correlation

cor(batch_patient_cluster_gzma$n, batch_patient_cluster_gzma$total, method = 'spearman')

ggplot(batch_patient_cluster_gzma, aes(x = total, y = n)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm') +
  theme_bw() +
  xlab('Total number of cells in patient') +
  ylab('Number of cells from cluster GZMA') +
  ggpubr::stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
                   label.x.npc = "left", label.y.npc = "top", size = 5)

cowplot::ggsave2('cluster_gzma_vs_total_correlation.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 300)

# Calculate frequencies between KRAS and mut

batch_patient_cluster_gzma_kras <- batch_patient_cluster_gzma %>% filter(KRAS != 'Unknown')

batch_patient_cluster_gzma_kras$freq_total <- batch_patient_cluster_gzma_kras$n / batch_patient_cluster_gzma_kras$total
batch_patient_cluster_gzma_kras$freq_treg <- batch_patient_cluster_gzma_kras$n / batch_patient_cluster_gzma_kras$total_treg

ggplot(batch_patient_cluster_gzma_kras, aes(x = KRAS, y = freq_total, fill = KRAS)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = 'jitter', size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  ylab('Frequency') +
  scale_fill_manual(values = col_vector[c(19,30,40)]) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             size = 8,
                             label = "p.signif", 
                             label.y = max(batch_patient_cluster_gzma_kras$freq_total, na.rm = TRUE) * 0.95,
                             label.x = 1.3) + NoLegend()

cowplot::ggsave2('S6D_cluster_type1_treg_freq_total.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 600)

ggplot(batch_patient_cluster_gzma_kras, aes(x = KRAS, y = freq_treg, fill = KRAS)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = 'jitter', size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  ylab('Frequency') +
  scale_fill_manual(values = col_vector[c(19,30,40)]) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             size = 8,
                             label = "p.signif", 
                             label.y = max(batch_patient_cluster_gzma_kras$freq_treg, na.rm = TRUE) * 0.95,
                             label.x = 1.3)

cowplot::ggsave2('3C_cluster_gzma_freq_tregs.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 600)

# Potential figure STUDY IFN-INDUCED CLUSTER------------------------------------------------------

cluster_of_interest <- 'PUT CLUSTER HERE'

batch_patient_stats <- integrated_tregs@meta.data %>% 
  group_by(seurat_clusters, orig.ident, Patient, KRAS) %>%
  summarize(n = n())

batch_patient_cluster_ifn <- batch_patient_stats %>% filter(seurat_clusters == cluster_of_interest) %>% filter(n >= 10)

patients_cluster_ifn <- batch_patient_cluster_ifn %>% pull(Patient)

patients_cluster_ifn_total_cells <- integrated@meta.data %>% filter(Patient %in% patients_cluster_ifn) %>% group_by(Patient) %>% summarize(total = n())

patients_cluster_ifn_total_tregs <- integrated_tregs@meta.data %>% filter(Patient %in% patients_cluster_ifn) %>% group_by(Patient) %>% summarize(total_treg = n())

batch_patient_cluster_ifn <- batch_patient_cluster_ifn %>%
  left_join(patients_cluster_ifn_total_cells, by = 'Patient') %>%
  left_join(patients_cluster_ifn_total_tregs, by = 'Patient')

# Check if the amount of this specific cells depend on the amount of total cells in Patient. In my point of view, this means that is can't find this cells 
# if the amount of total cells in Patient is low if there is an correlation

cor(batch_patient_cluster_ifn$n, batch_patient_cluster_ifn$total, method = 'spearman')

ggplot(batch_patient_cluster_ifn, aes(x = total, y = n)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm') +
  theme_bw() +
  xlab('Total number of cells in patient') +
  ylab('Number of cells from cluster ifn') +
  ggpubr::stat_cor(method = "pearson", aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
                   label.x.npc = "left", label.y.npc = "top", size = 5)

cowplot::ggsave2('cluster_ifn_vs_total_correlation.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 300)

# Calculate frequencies between KRAS and mut

batch_patient_cluster_ifn_kras <- batch_patient_cluster_ifn %>% filter(KRAS != 'Unknown')

batch_patient_cluster_ifn_kras$freq_total <- batch_patient_cluster_ifn_kras$n / batch_patient_cluster_ifn_kras$total
batch_patient_cluster_ifn_kras$freq_treg <- batch_patient_cluster_ifn_kras$n / batch_patient_cluster_ifn_kras$total_treg

ggplot(batch_patient_cluster_ifn_kras, aes(x = KRAS, y = freq_total, fill = KRAS)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = 'jitter', size = 3) +
  theme_bw() +
  scale_fill_manual(values = col_vector[c(19,30,40)])

cowplot::ggsave2('cluster_ifn_freq_total.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 300)

ggplot(batch_patient_cluster_ifn_kras, aes(x = KRAS, y = freq_treg, fill = KRAS)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = 'jitter', size = 3) +
  theme_bw() +
  scale_fill_manual(values = col_vector[c(19,30,40)])

cowplot::ggsave2('cluster_ifn_freq_tregs.pdf', path = figures_folder, width = 15, height = 10, units = 'cm', dpi = 300)

# Potential figure STUDY TREG GZMA CLONALITY------------------------------------------------------

integrated_tregs_th1lile <- subset(integrated_tregs, subset = seurat_clusters == 10)

intgr <- integrated_tregs_th1lile

intgr@meta.data$patient_id <- paste0(intgr@meta.data$orig.ident, '_', intgr@meta.data$Patient)
intgr@meta.data$Subset <- intgr@meta.data$seurat_clusters
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$Barcode <- row.names(intgr[[]])

#  No less than N cells in each Subset in each patient

N = 20
intgr_filtered_cells <- group_by(intgr[[]], Subset, patient_id) %>% filter(n() > N) 
intgr <- subset(intgr, cells = intgr_filtered_cells$Barcode)

intgr$Barcode <- NULL
intgr$barcode <- NULL # in case there is such column

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  tibble::rownames_to_column(var = "barcode")

umap_tx %>%
  dplyr::select(patient_id, Subset, TCR2nt, KRAS) %>%
  tidyr::drop_na() %>% # we take only TRB containing cells %>%
  group_by(patient_id, Subset, KRAS) %>%
  mutate(n_cells_in_cluster = n()) %>%
  distinct() %>% # keep only unique clonotypes
  group_by(patient_id, Subset, KRAS, n_cells_in_cluster) %>%
  summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
  group_by(patient_id, KRAS) %>%
  mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
         fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
         n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

clonality_data %>%
  mutate(ratio = n_clones_in_cluster / n_cells_in_cluster) -> clonality_data

clonality_data$Subset <- 'T reg GZMA+ CCL5+'

plotClonalityCelltypes <- function(x){
  # Shows clonality between cell types without spliiting into groups by covariate
  clonality_data <- x
  c1 <- ggplot(clonality_data, aes(x = Subset, y = ratio, fill = Subset))  +
    coord_flip() +
    geom_point(key_glyph = "point",position=position_jitterdodge(), size = 5) +
    theme_minimal() +
    ggtitle("TCR\u03b2 clonality") +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
    #scale_color_identity(guide = "legend", labels = names(subsets_cols), breaks = subsets_cols) +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          panel.grid.major.y = element_blank()) +
    #scale_x_log10() +
    #scale_y_log10() +
    xlab("") +
    ylab("N TCR\u03b2 clonotypes in cluster / N cells in cluster")  + 
    NoLegend() + 
    ylim(0, 1)
  
  c2 <- ggplot(data = clonality_data, aes(x = Subset)) +
    coord_flip() +
    geom_bar(stat = 'count', position = 'dodge', fill = '#99A3A4') +
    xlab('') +
    ylab('') +
    theme_bw() +
    NoLegend() + 
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  fig <- ggpubr::ggarrange(c1, c2, nrow = 1, align = 'h', widths = c(5, 1), vjust = 0.1)
  
  plot(fig)
}

plot_clonality <- plotClonalityCelltypes(clonality_data)

ggsave(plot = plot_clonality, paste0(figures_folder, "tregs_gzma_clonality_20cells.pdf"), bg="white", height = 4, width = 7, dpi = 300)

# 4A Clonality -------------------------------------------------------------------


library(Seurat)
library(tidyverse)
library(here)
library(pals)
library(ggpubr)
library(RColorBrewer)

setwd('/projects/sle_jul_23_gabibov/')


source('/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/luad_figures_tables/colors_fixed.R')
source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')

subsets_cols <- subsets_cols_level2

# Annotations should be stored under metadata Subset column, TCR cdr3 beta under TCR2nt,
# donor under patient_id

figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed"
data_folder <- '/projects/sle_jul_23_gabibov/luad/tcr_analysis/data/'

path_to_intgr_seurat <- "/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds"

intgr_full <- read_rds(path_to_intgr_seurat)


integrated_backup <- intgr_full

intgr_full@meta.data$CelltypeAnnotationLow <- intgr_full@meta.data$CelltypeAnnotationLevel2


intgr_full <- subset(intgr_full, subset = CelltypeAnnotationLevel1  != 'Unknown')


intgr <- intgr_full
#intgr <- subset(intgr_full, subset = KRAS != 'Unknown')


intgr@meta.data$patient_id <- intgr@meta.data$Patient
intgr@meta.data$Subset <- intgr@meta.data$CelltypeAnnotationLow
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$Barcode <- row.names(intgr[[]])

#  Delete cells without TCR data
cells_keep <- rownames(intgr@meta.data)[!is.na(intgr@meta.data$TCR2nt)]
intgr <- subset(intgr, cells = cells_keep)
#  No less than N cells in each Subset in each patient

N = 50
intgr_filtered_cells <- group_by(intgr[[]], Subset, patient_id) %>% filter(n() > N) 
intgr <- subset(intgr, cells = intgr_filtered_cells$Barcode)

intgr$Barcode <- NULL
intgr$barcode <- NULL # in case there is such column

# extract UMAP coordinates for cells 

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode")

# clonality plot

# 100% - all TRB clones

# umap_tx %>%
#   dplyr::select(patient_id, Subset, TCR2nt) %>%
#   drop_na() %>% # we take only TRB containing cells %>%
#   group_by(patient_id, Subset) %>%
#   mutate(n_cells_in_cluster = n()) %>%
#   distinct() %>% # keep only unique clonotypes
#   group_by(patient_id, Subset,  n_cells_in_cluster) %>%
#   summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
#   group_by(patient_id) %>%
#   mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
#          fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
#          n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

umap_tx %>%
  dplyr::select(patient_id, Subset, TCR2nt, KRAS) %>%
  drop_na() %>% # we take only TRB containing cells %>%
  group_by(patient_id, Subset, KRAS) %>%
  mutate(n_cells_in_cluster = n()) %>%
  distinct() %>% # keep only unique clonotypes
  group_by(patient_id, Subset, KRAS, n_cells_in_cluster) %>%
  summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
  group_by(patient_id, KRAS) %>%
  mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
         fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
         n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

clonality_data %>%
  mutate(ratio = n_clones_in_cluster / n_cells_in_cluster) %>%
  mutate(Subset = factor(Subset, ordered = T,
                         levels=names(subsets_cols))) -> clonality_data

clonality_data$colors <- subsets_cols[clonality_data$Subset]

# # Filter out cell types with less than 3 patients UNIQUE FEATURE
# 
# patients_per_celltypes <- clonality_data %>%
#   distinct(patient_id, KRAS, Subset) %>%
#   group_by(Subset, KRAS) %>%
#   summarise(n_patients = n(), .groups = "drop") %>%
#   ungroup %>%
#   group_by(Subset) %>%
#   filter(all(n_patients[KRAS == "WT"] >= 3),
#          all(n_patients[KRAS == "mut"] >= 3),
#          all(n_patients[KRAS == "Unknown"] >= 3)) %>%
#   ungroup()
# 
# clonality_data <- clonality_data %>%
#   filter(Subset %in% unique(patients_per_celltypes$Subset))
# 
# # Filter out cell types with less than 3 patients UNIQUE FEATURE

# # FILTER BY MEDIAN ABSOLUTE DEVIATION
# 
# stats <- clonality_data %>%
#   group_by(Subset, KRAS) %>%
#   summarise(
#     med = median(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones),
#     mad_val = mad(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones),
#     .groups = "drop"
#   )
# 
# clonality_data <- clonality_data %>%
#   left_join(stats, by = c("Subset", "KRAS")) %>%
#   mutate(
#     mad_score = abs(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones - med) / mad_val
#   )
# 
# # PLOT DISTRIBUTION OF MAD SCORES
# 
# threshold <- quantile(clonality_data$mad_score, 0.9, na.rm = TRUE)
# 
# mad_score_distribution <- ggplot(clonality_data, aes(x = mad_score)) +
#   geom_histogram(bins = 105) +
#   theme_bw() +
#   xlab('MAD score') +
#   ylab('Count') +
#   geom_vline(xintercept = threshold, color = 'red', linetype = 'dashed')
# 
# mad_score_distribution
# 
# ggsave(plot = mad_score_distribution, paste0(figures_folder, "madScoreFull_prepostT_09perc_min50cells.png"), bg="white", height = 8, width = 10, dpi = 300, units = "cm")
# 
# 
# # FILTER OUT ONLY SAMPLES WITH HIGH MAD
# 
# clonality_data <- clonality_data %>%
#   filter(mad_score <= threshold)

# FILTER BY MEDIAN ABSOLUTE DEVIATION

clonality_data$clonality <- 1 - clonality_data$ratio 

# write.csv(clonality_data, paste0(data_folder, 'clonalityFullNSCLC_Level2.csv'), row.names = FALSE)

plotClonalityCelltypes <- function(x){
  # Shows clonality between cell types without spliiting into groups by covariate
  clonality_data <- x
  c1 <- ggplot(clonality_data, aes(x = Subset, y = clonality, color = colors))  +
    coord_flip() +
    # geom_boxplot(alpha = 0.7) +
    # stat_compare_means(method = "anova", show.legend = FALSE,
    #                    label.x = length(unique(clonality_data$Subset)),
    #                    label.y = 0.5,
    #                    size = 3) +
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = paste0(clonality_data %>% group_by(Subset) %>% summarise(mean = mean(ratio)) %>% arrange(mean) %>% tail(1) %>% pull(Subset)),
    #                    label.y = 0, hide.ns = T, p.adjust.method = "BH") +
    geom_point(key_glyph = "point",position=position_jitterdodge(), size = 5) +
    theme_minimal() +
    ggtitle("TCR\u03b2 clonality") +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
    scale_color_identity(guide = "legend", labels = names(subsets_cols), breaks = subsets_cols) +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          panel.grid.major.y = element_blank()) +
    #scale_x_log10() +
    #scale_y_log10() +
    xlab("") +
    ylab("Clonality")  + 
    NoLegend() + 
    ylim(0, 1.0) +
    scale_y_reverse(breaks = seq(0, 1, 0.25))
  
  c2 <- ggplot(data = clonality_data, aes(x = Subset)) +
    coord_flip() +
    geom_bar(stat = 'count', position = 'dodge', fill = '#99A3A4') +
    xlab('') +
    ylab('') +
    theme_bw() +
    NoLegend() + 
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_blank())
  
  fig <- ggpubr::ggarrange(c1, c2, nrow = 1, align = 'h', widths = c(5, 1), vjust = 0.1)
  
  plot(fig)
}

plot_clonality <- plotClonalityCelltypes(clonality_data)

plot_clonality

ggsave(plot = plot_clonality, paste0(figures_folder, "/4A_clonality_full_prepostTreatments_NoMadF_min50cells_Level2.pdf"), bg="white", height = 7, width = 10, dpi = 600)

gc()

# 4B blood_nsclc -------------------------------------------------------------



library(Seurat)
library(tidyverse)
library(here)
library(pals)
library(ggpubr)
library(RColorBrewer)

setwd('/projects/sle_jul_23_gabibov/')


source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

# subsets_cols <- subsets_cols_level2

# Annotations should be stored under metadata Subset column, TCR cdr3 beta under TCR2nt,
# donor under patient_id

data_folder <- '/projects/sle_jul_23_gabibov/luad/tcr_analysis/data/'

clonality_blood <- read.csv(paste0(data_folder, 'clonalityFullBlood.csv'))
clonality_blood$Status <- 'Blood'

clonality_nsclc <- read.csv(paste0(data_folder, 'clonalityFullNSCLC_Level1.csv'))
clonality_nsclc$KRAS <- NULL
clonality_nsclc$Status <- 'NSCLC'

#-----------------------------------------------------------------RENAME SUBSETS IN NSCLC TO BLOOD

subset_rename_map <- c(
  "Naive"                = "Tcm/Naive",
  "Naive_RTE"            = "Tcm/Naive",
  "CentMem1"             = "Tcm/Naive",
  "CentMem2"             = "Tcm/Naive",
  "EffMem_Th17"          = "Tem Th17 like",
  "EffMem_Th1"           = "Tem Th1 like",
  "Treg"                 = "T reg",
  "Cycling"              = "Cycling",
  "PD1high"              = "PD1 high",
  "Temra_cytotoxic_Th1"  = "Cytotoxic TEMRA",
  "Tfh"                  = "Tem Tfh like",
  "EffMem_IFN_induced"   = "IFN-induced"
)

clonality_blood <- clonality_blood %>%
  mutate(Subset = dplyr::recode(Subset, !!!subset_rename_map))

# subset_rename_map <- c(
#   "Tcm/Naive"            = "Naive",
#   "Tem Th17 like"        = "EffMem_Th17",
#   "Tem Th1 like"         = "EffMem_Th1",
#   "Trm Th1 like"         = "EffMem_Th1",
#   "T reg"                = "Treg",
#   "Cycling"              = "Cycling",
#   "PD1 high"             = "PD1high",
#   "Cytotoxic PD1 high"   = "PD1high",
#   "Cytotoxic TEMRA"      = "Temra_cytotoxic_Th1",
#   "Tem Tfh like"         = "Tfh",
#   "Trm Tfh like"         = "Tfh",
#   "IFN-induced"          = "IFN_induced"
# )
# 
# clonality_nsclc <- clonality_nsclc %>%
#   mutate(Subset = recode(Subset, !!!subset_rename_map))

clonality_full <- rbind(clonality_blood, clonality_nsclc)

clonality_full <- clonality_full %>% filter(Subset %in% intersect(unique(clonality_nsclc$Subset), unique(clonality_blood$Subset)))

clonality_full$Subset <- factor(clonality_full$Subset, levels = c('Tcm/Naive',
                                                                  'Tem Th17 like',
                                                                  'Tem Th1 like',
                                                                  'Cytotoxic TEMRA',
                                                                  'Tem Tfh like',
                                                                  'PD1 high', 
                                                                  'IFN-induced',
                                                                  'T reg'))

celltypes_to_delete <- clonality_full %>% group_by(Subset, Status) %>% summarise(count = n(), median_count = median(clonality)) %>% filter(count < 4) %>% pull(Subset)

clonality_full <- clonality_full %>% filter(!Subset %in% celltypes_to_delete)

#-----------------------------------------------------------------STATISTICAL TESTING

wilcox_results <- clonality_full %>%
  mutate(Status = as.factor(Status)) %>%
  group_by(Subset) %>%
  summarise(
    p = tryCatch(
      wilcox.test(clonality ~ Status, exact=FALSE)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  mutate(star_adj = ifelse(p_adj < 0.05, "*", ""))

clonality_full <- clonality_full %>%
  left_join(wilcox_results, by = "Subset")


library(car)

levene_results <- clonality_full %>%
  mutate(Status = as.factor(Status)) %>%
  group_by(Subset) %>%
  summarise(
    p = tryCatch(
      leveneTest(clonality ~ Status)$`Pr(>F)`[1],
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    star_adj = ifelse(p_adj < 0.05, "*", "")
  )


#-----------------------------------------------------------------PLOTTING

plotClonalityGroups <- function(x){
  # Shows clonality between groups inside covariate
  clonality_data <- x
  c1 <- ggplot(clonality_data, aes(x = Subset, y = clonality, fill = Status))  +
    coord_flip() +
    geom_point(key_glyph = "point",position=position_jitterdodge(), size = 1, alpha = 0.4) +
    geom_boxplot(alpha = 1, outliers = FALSE) +
    # stat_compare_means(method = "anova", show.legend = FALSE, 
    #                    label.x = length(unique(clonality_data$Subset)), 
    #                    label.y = 0.5,
    #                    size = 3) +
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = paste0(max(clonality_data$Subset)),
    #                    label.y = 0, hide.ns = T) +
    theme_minimal() +
    ggtitle("TCR\u03b2 clonality") +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
    scale_color_identity(guide = "legend") +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          panel.grid.major.y = element_blank()) +
    #scale_x_log10() +
    #scale_y_log10() +
    xlab("") +
    ylab("Clonality")  + 
    NoLegend() + 
    ylim(0, 1) +
    scale_fill_manual(values = col_vector[c(25,30)]) +
    scale_y_reverse(breaks = seq(0, 1, 0.25)) +
    geom_text(
      data = wilcox_results,
      aes(y = 1, x = Subset, label = star_adj),
      inherit.aes = FALSE,
      vjust = 0.75,
      size = 8,
      color = 'darkgreen'
    )
  
  c2 <- ggplot(data = clonality_data, aes(x = Subset, fill = Status)) +
    coord_flip() + #ylim = c(0, 20)
    geom_bar(stat = 'count', position = 'dodge') +
    xlab('') +
    ylab('') +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.y = element_blank()) +
    #scale_y_continuous(breaks = seq(0, 20, by = 5)) +
    scale_fill_manual(values = col_vector[c(25, 30)])
  
  fig <- ggpubr::ggarrange(c1, c2, nrow = 1, align = 'h', common.legend = TRUE, legend = 'bottom', widths = c(5, 1), vjust = 0.1)
  
  plot(fig)
}

plot_clonality <- plotClonalityGroups(clonality_full)

plot_clonality

ggsave(plot = plot_clonality, paste0(figures_folder, "/4B_clonality_full_prepostTreatments_NoMadF_min50cells_Blood_NSCLC_2.pdf"), bg="white", height = 7, width = 10, dpi = 600)


# 4D D_metric_merged ---------------------------------------------------------


library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(forcats)
library(RColorBrewer)

# Annotations should be stored under metadata Subset column, TCR cdr3 beta under TCR2nt,
# donor under patient_id

setwd('/projects/sle_jul_23_gabibov/')


source('/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/luad_figures_tables/colors_fixed.R')
source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')

subsets_cols <- subsets_cols_level1

N <- 50
K <- 3

figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed"
data_folder <- '/projects/sle_jul_23_gabibov/luad/tcr_analysis/data/'

mat_mut <- readRDS(paste0(data_folder, paste0('matrix_krasMut_prepostTreatment_', N, 'clones_', K, 'patientPerCelltype_Level1_PD1merged.rds')))
mat_wt <- readRDS(paste0(data_folder, paste0('matrix_krasWT_prepostTreatment_', N, 'clones_', K, 'patientPerCelltype_Level1_PD1merged.rds')))
mat_full <- readRDS(paste0(data_folder, paste0('matrix_krasFull_prepostTreatment_', N, 'clones_', K, 'patientPerCelltype_Level1_PD1merged.rds')))

#-----------------------------------------------------------------COMBINE MATRICES

common_celltypes <- intersect(rownames(mat_mut), rownames(mat_wt))
common_celltypes <- intersect(common_celltypes, rownames(mat_full))

mat_mut_sub <- mat_mut[common_celltypes, common_celltypes]
mat_wt_sub  <- mat_wt[common_celltypes, common_celltypes]
mat_full_sub <- mat_full[common_celltypes, common_celltypes]


# Combine 2 matrix
# combined_mat <- mat_mut_sub
# combined_mat[lower.tri(combined_mat)] <- mat_wt_sub[lower.tri(mat_wt_sub)]

#mat_plot <- mat_wt_sub

mat_plot <- mat_mut_sub
mat_plot[upper.tri(mat_plot)] <- NA

diag(mat_plot) <- -1

max_value <- 0.0004

# # Get all unique cell types
# all_celltypes <- union(rownames(mat_mut), rownames(mat_wt))
# 
# # Initialize full-size matrices with NA
# mat_mut_full <- matrix(NA, nrow = length(all_celltypes), ncol = length(all_celltypes),
#                        dimnames = list(all_celltypes, all_celltypes))
# mat_wt_full <- mat_mut_full  # same structure
# 
# # Fill available values
# mat_mut_full[rownames(mat_mut), colnames(mat_mut)] <- mat_mut
# mat_wt_full[rownames(mat_wt), colnames(mat_wt)] <- mat_wt
# 
# # Combine into one matrix
# combined_mat <- mat_mut_full  # start with one matrix
# combined_mat[lower.tri(combined_mat)] <- mat_wt_full[lower.tri(combined_mat)]

col_fun = colorRamp2(c(-1, 0, max_value), c("grey", "white", "red"))

clust_order <- names(subsets_cols)

cluster_cols <- subsets_cols

subsets <- colnames(mat_plot)[match(clust_order, colnames(mat_plot), nomatch = 0L)]

p1 <- Heatmap(mat_plot,
              na_col = 'white',
              column_title = expression(" "),
              column_order = subsets,
              column_names_side = "bottom",
              column_names_rot = 45,
              row_order = subsets,
              show_row_names = T,
              row_names_side = "left",
              width = unit(16.6, "cm"),
              height = unit(16.6, "cm"),
              #rect_gp = gpar(col = "gray", lwd = 1),
              col = col_fun,
              heatmap_legend_param = list(
                title = expression("log D-metric"),
                at = c(0, max_value/2, max_value),
                labels = c("0", paste0(max_value/2), paste0(max_value))
              ),
              bottom_annotation = HeatmapAnnotation(clusters = subsets,
                                                    col = list(clusters = cluster_cols),
                                                    show_annotation_name = F,
                                                    show_legend = F),
              left_annotation = rowAnnotation(clusters = subsets,
                                              col = list(clusters = cluster_cols),
                                              show_legend = F,
                                              show_annotation_name = F))
plot(p1)

pdf(file=paste0(figures_folder, paste0("/4D_D_metricTriangle_mergedMut_prepostTreatment_", N, "clonesF_", K, "patientPerCelltype_Level1_PD1high.pdf")), width=11,height=11)
draw(p1)
dev.off()
pdf(file=paste0(figures_folder, paste0("/4D_D_metricTriangle_mergedWT_prepostTreatment_", N, "clonesF_", K, "patientPerCelltype_Level1_PD1high.pdf")), width=11,height=11)
draw(p1)
dev.off()
rm(mat_plot)


# 4С Chord -------------------------------------------------------------------


library(Seurat)
library(tidyverse)
library(here)
library(pals)
library(ggpubr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(forcats)

setwd('/projects/sle_jul_23_gabibov/')


source('/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/luad_figures_tables/colors_fixed.R')
source('/projects/sle_jul_23_gabibov/luad/helpers/functions.R')

subsets_cols <- subsets_cols_level1

# Annotations should be stored under metadata Subset column, TCR cdr3 beta under TCR2nt,
# donor under patient_id

figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras/colors_fixed"
data_folder <- '/projects/sle_jul_23_gabibov/luad/tcr_analysis/data/'

path_to_intgr_seurat <- "/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds"

intgr_full <- read_rds(path_to_intgr_seurat)

#-----------------------------------------------------------------ASSIGN CORRECT VERSION OF ANNOTATION

integrated_backup <- intgr_full

intgr_full@meta.data$CelltypeAnnotationLow <- intgr_full@meta.data$CelltypeAnnotationLevel1

#-----------------------------------------------------------------DELETE UNKNOWN AND ARTEFACT CELL TYPES

intgr_full <- subset(intgr_full, subset = CelltypeAnnotationLevel1 != 'Unknown')

#-----------------------------------------------------------------SUBSET BY KRAS STATUS

intgr <- intgr_full

#-----------------------------------------------------------------CHANGE SPECIFIC CELL TYPES

intgr@meta.data$CelltypeAnnotationLow <- dplyr::recode(intgr@meta.data$CelltypeAnnotationLow, 'Cytotoxic PD1 high' = 'PD1 high')

#-----------------------------------------------------------------CREATE REQUIRED COLUMNS

intgr@meta.data$patient_id <- intgr@meta.data$Patient
intgr@meta.data$Subset <- intgr@meta.data$CelltypeAnnotationLow
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$TCR2nt[intgr$TCR2nt == "NA"] <- NA
tcr_obj_cells <- dplyr::select(intgr[[]], TCR2nt) %>% drop_na()
intgr <- subset(intgr, cells = row.names(tcr_obj_cells))

intgr$donor_TCR2nt <- paste0(intgr$patient_id, "_", intgr$TCR2nt) # Count only intra-patient plasticity


meta <- intgr@meta.data

clone_celltype <- meta %>%
  drop_na(cdr3_nt2) %>%
  distinct(cdr3_nt2, CelltypeAnnotationLow)

# Find all combinations of cell types sharing the same clone. Delete combinations for the same cell types
clone_pairs <- inner_join(clone_celltype, clone_celltype, by = "cdr3_nt2") #%>%
#filter(CelltypeAnnotationLow.x != CelltypeAnnotationLow.y)

# Calculate number of shared clones
clone_links <- clone_pairs %>%
  dplyr::count(CelltypeAnnotationLow.x, CelltypeAnnotationLow.y, name = "shared_clones")

# Account for symmetry
clone_links_sym <- clone_links %>%
  rowwise() %>%
  mutate(
    pair = list(sort(c(CelltypeAnnotationLow.x, CelltypeAnnotationLow.y)))
  ) %>%
  ungroup() %>%
  mutate(
    from = sapply(pair, `[`, 1),
    to   = sapply(pair, `[`, 2)
  ) %>%
  group_by(from, to) %>%
  summarise(shared_clones = sum(shared_clones), .groups = "drop")

# Add value for inter-cluster links
clone_links <- clone_links %>%
  mutate(intercluster = case_when(
    CelltypeAnnotationLow.x == CelltypeAnnotationLow.y ~ 0,
    TRUE ~ 1
  ))

clone_links_sym <- clone_links_sym %>%
  mutate(intercluster = case_when(
    from == to ~ 0,
    TRUE ~ 1
  ))


pdf(paste0(figures_folder, "/4C_chordPlot2_Full_scaled_Level1.pdf"))

# Your existing chordDiagram with scaled values
chordDiagram(clone_links,
             grid.col = subsets_cols,
             directional = 1,
             annotationTrack = c("grid"),
             scale = TRUE,
             preAllocateTracks = list(track.height = mm_h(5)),
             annotationTrackHeight = mm_h(5))

# First axis: scaled percentages (0%–100%)
for(si in get.all.sector.index()) {
  circos.axis(h = "top", 
              sector.index = si, 
              major.at = seq(0, 1, by = 0.25),
              labels.cex = 0.35,
              minor.ticks = 4,
              major.tick.length = convert_y(1.5, "mm"))
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), -0.75, si, sector.index = si, track.index = 1, 
              facing = "bending.inside", col = 'black', niceFacing = TRUE, cex = 0.6)
}

# Prepare raw counts per sector
clone_vector <- clone_celltype %>% 
  group_by(CelltypeAnnotationLow) %>% 
  summarise(count = n()) %>% 
  deframe()

# Second axis (also above sectors): unscaled (raw clone counts)
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  max_val = clone_vector[sector.name]
  
  # Add dotted reference lines above sectors
  circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3, col = "grey50")
  
  # Define positions for raw count axis ticks
  ticks_at = seq(0, 1, by = 0.25) 
  labels = round(seq(0, max_val, length.out = length(ticks_at)))
  
  for(i in seq_along(ticks_at)) {
    if (labels[i] != 0) { 
      circos.text(xlim[1] + diff(xlim) * ticks_at[i],
                  mean(ylim) + 0.3,
                  labels[i], 
                  cex = 0.3, adj = c(0.2, 0.5), facing = 'clockwise', niceFacing = TRUE)
    }
  }
  
}, bg.border = NA)

circos.clear()

dev.off()



# S5B ---------------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(RColorBrewer)
library(kableExtra)
library(cowplot)
library(miloR)
library(SingleCellExperiment)

setwd('/projects/sle_jul_23_gabibov/')


source('luad/helpers/colors.R')
source('luad/helpers/functions.R')


print('Read data')

figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras"
data_folder <- '/projects/sle_jul_23_gabibov/luad/milo/data/'

integrated <- readRDS('/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds')

#-----------------------------------------------------------------DELETE UNKNOWN AND ARTEFACT CELL TYPES-

integrated <- subset(integrated, subset = CelltypeAnnotationLevel1 != 'Unknown')

#-----------------------------------------------------------------FILTER BY KRAS STATUS

#integrated_kras <- subset(integrated, subset = KRAS != 'Unknown')

integrated_kras <- integrated

#-----------------------------------------------------------------CREATE MILO OBJECT

print('Create Milo object')
#integrated[["RNA"]] <- as(integrated[["RNA"]], Class="Assay")

integrated_kras_sce <- as.SingleCellExperiment(integrated_kras, assay = 'RNA')

integrated_kras_milo <- Milo(integrated_kras_sce)

#-----------------------------------------------------------------CONSTRUCT KNN GRAPH

print('Construct KNN graph')

integrated_kras_milo <- buildGraph(integrated_kras_milo, k = 30, d = 30, reduced.dim = 'INTEGRATED.RPCA')

#-----------------------------------------------------------------NHOODS

print('Calculate Nhoods')

integrated_kras_milo <- makeNhoods(integrated_kras_milo, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = 'INTEGRATED.RPCA')

plotNhoodSizeHist(integrated_kras_milo)

cowplot::ggsave2('S5B_nhoods_size.pdf', path = figures_folder, height = 15, width = 20, units = 'cm', dpi = 300)

#-----------------------------------------------------------------COUNT CELLS IN NBHOODS

print('Count cells')

integrated_kras_milo <- countCells(integrated_kras_milo, meta.data = as.data.frame(colData(integrated_kras_milo)), sample="Patient")

#----------------------------------------------------------------DESIGN MATRIX

print('Create matrix')

design <- data.frame(colData(integrated_kras_milo))[,c("Patient", "KRAS", "orig.ident")]

## Convert batch info from integer to factor
design$orig.ident <- as.factor(design$orig.ident) 
design <- distinct(design)
rownames(design) <- design$Patient

design

#----------------------------------------------------------------NBHOOD CONNECTIVITY

print('Calculate connectivity')

integrated_kras_milo <- calcNhoodDistance(integrated_kras_milo, d=30, reduced.dim = "INTEGRATED.RPCA")

#----------------------------------------------------------------DIFFERENTIAL ABUNDANCE TESTING

print('Testing')

contrast <- c("KRASmut - KRASWT")

da_results <- testNhoods(integrated_kras_milo, 
                         design = ~ 0 + KRAS,
                         design.df = design, 
                         model.contrasts = contrast,
                         fdr.weighting="graph-overlap", 
                         norm.method="TMM")

da_results_batch <- testNhoods(integrated_kras_milo, 
                               design = ~ orig.ident + 0 + KRAS,
                               design.df = design, 
                               model.contrasts = contrast,
                               fdr.weighting="graph-overlap", 
                               norm.method="TMM")

da_results_patient <- testNhoods(integrated_kras_milo, 
                                 design = ~ Patient + 0 + KRAS,
                                 design.df = design, 
                                 model.contrasts = contrast,
                                 fdr.weighting="graph-overlap", 
                                 norm.method="TMM")

#----------------------------------------------------------------SAVE MILO

print('Save results')

#saveRDS(integrated_kras_milo, paste0(data_folder, 'integrated_milo_final.rds'))
#write.csv(da_results, paste0(data_folder, 'da_milo_final.csv'))
#write.csv(da_results_batch, paste0(data_folder, 'da_milo_final_batch.csv'))
#write.csv(da_results_patient, paste0(data_folder, 'da_milo_final_patient.csv'))

#----------------------------------------------------------------READ MILO

integrated_kras_milo <- readRDS(paste0(data_folder, 'integrated_milo_final.rds'))
da_results <- read.csv(paste0(data_folder, 'da_milo_annotated.csv'))

#----------------------------------------------------------------ASSIGN NHOODS TO DISCRETE CLUSTERS

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)


da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel1")
da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel2")

da_results$CelltypeAnnotationLevel1 <- factor(da_results$CelltypeAnnotationLevel1, levels = names(subsets_cols_level1))
da_results$CelltypeAnnotationLevel2 <- factor(da_results$CelltypeAnnotationLevel2, levels = setdiff(names(subsets_cols_level2), 'Unknown'))

plotDAbeeswarm(da_results, group.by = "CelltypeAnnotationLevel2")

cowplot::ggsave2('da_plot_final.pdf', path = figures_folder, height = 20, width = 20, units = 'cm', dpi = 300)

# 2D  FILTER NHOODS BY FREQUENCY IN CLUSTERS------------------------------------------------------------

nhood_freq <- 0.95

fdr_threshold <- 0.05

da_results_filtered <- filter(da_results, CelltypeAnnotationLevel2_fraction >= nhood_freq)

da_medians <- da_results_filtered %>%
  group_by(CelltypeAnnotationLevel2) %>%
  summarise(median_logFC = median(logFC, na.rm = TRUE)) %>%
  mutate(group_by = CelltypeAnnotationLevel2)

da_medians$group_by <- factor(da_medians$group_by, levels = unique(da_results_filtered$CelltypeAnnotationLevel2))

da_medians$pos_y <- seq_along(da_medians$group_by)

plotDAbeeswarm(da_results_filtered, group.by = "CelltypeAnnotationLevel2", alpha = fdr_threshold) +
  ggplot2::geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
  geom_point(data = da_medians,
             aes(y = median_logFC, x = pos_y),
             inherit.aes = FALSE,
             color = "red", size = 3) +
  geom_text(data = da_medians,
            aes(y = -4.1, x = pos_y, label = round(median_logFC, 2)),
            inherit.aes = FALSE,
            size = 5,
            hjust = 1) +
  scale_y_continuous(
    limits = c(-5, 4),
    breaks = c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
    minor_breaks = NULL
  ) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

cowplot::ggsave2('2D_da_plot_final_filtered.pdf', path = figures_folder, height = 20, width = 25, units = 'cm', dpi = 300)

#----------------------------------------------------------------GOODNESS-OF-FIT TESTS

# The idea here is to test that distribution of number of statistically significant nhoods in KRAS-mut (positive) and KRAS-WT (negative) differs from equal 0.5:0.5

metrics_summary <- da_results_filtered %>%
  mutate(
    is_signif = SpatialFDR < fdr_threshold,
    logFC_dir = case_when(
      is_signif & logFC > 0 ~ "signif_pos",
      is_signif & logFC < 0 ~ "signif_neg",
      TRUE ~ "other"
    )
  ) %>%
  group_by(CelltypeAnnotationLevel2) %>%
  summarise(
    total_nhoods = n(),
    n_significant = sum(is_signif),
    freq_significant = mean(is_signif),
    n_nonsignificant = sum(!is_signif),
    freq_nonsignificant = mean(!is_signif),
    n_significant_logFC_pos = sum(logFC_dir == "signif_pos"),
    n_significant_logFC_neg = sum(logFC_dir == "signif_neg"),
    freq_logFC_pos_total = n_significant_logFC_pos / total_nhoods,
    freq_logFC_neg_total = n_significant_logFC_neg / total_nhoods,
    freq_logFC_pos_within_signif = n_significant_logFC_pos / n_significant,
    freq_logFC_neg_within_signif = n_significant_logFC_neg / n_significant,
    
    # Chi-square goodness-of-fit test
    chisq_gof_p = purrr::map2_dbl(n_significant_logFC_pos, n_significant_logFC_neg, ~ {
      obs <- c(.x, .y)
      if (sum(obs) > 0) {
        chisq.test(obs, p = c(0.5, 0.5))$p.value
      } else {
        NA_real_
      }
    }),
    
    # Fisher's exact test version of the same logic (observed vs expected counts)
    fisher_gof_p = purrr::map2_dbl(n_significant_logFC_pos, n_significant_logFC_neg, ~ {
      total <- .x + .y
      if (total > 0) {
        expected <- total / 2
        tbl <- matrix(c(.x, .y, expected, expected), nrow = 2)
        fisher.test(tbl)$p.value
      } else {
        NA_real_
      }
    })
  ) %>%
  # Apply multiple testing correction
  mutate(
    chisq_gof_p_adj_bh = p.adjust(chisq_gof_p, method = "BH"),
    fisher_gof_p_adj_bh = p.adjust(fisher_gof_p, method = "BH"),
    chisq_gof_p_adj_bonf = p.adjust(chisq_gof_p, method = "bonferroni"),
    fisher_gof_p_adj_bonf = p.adjust(fisher_gof_p, method = "bonferroni")
  )

#write.csv(metrics_summary, paste0(data_folder, 'metrics_summary_fdr005_freq09.csv'))

# PseudoS5C--PLOT NHOODS ON KNN GRAPH------------------------------------------------------------

integrated_kras_milo <- buildNhoodGraph(integrated_kras_milo)

combined_plot <- DimPlot(integrated_kras, group.by = 'CelltypeAnnotationLevel2', raster = FALSE, cols = subsets_cols_level2, pt.size = 2) + 
  plotNhoodGraphDA(integrated_kras_milo, da_results, alpha=0.05, size_range=c(2,8), highlight.da = 0) +
  plot_layout(guides = "collect")

cowplot::ggsave2(filename = file.path(figures_folder, 'umaps_final.png'),
                 plot = combined_plot,
                 height = 25, width = 50, units = 'cm', dpi = 300)

sig_nhoods <- as.character(da_results[da_results$SpatialFDR < fdr_threshold, "Nhood"])

combined_plot2 <- DimPlot(integrated_kras, group.by = 'CelltypeAnnotationLevel2', raster = FALSE, cols = subsets_cols_level2) + 
  plotNhoodGraphDA(integrated_kras_milo, da_results, subset.nhoods = sig_nhoods, alpha=0.05) +
  plot_layout(guides = "collect")

cowplot::ggsave2(filename = file.path(figures_folder, 'pseudoS5C_umaps_final_filtrered.pdf'),
                 plot = combined_plot2,
                 height = 30, width = 60, units = 'cm', dpi = 300)


library(Seurat)
library(dplyr)
library(DescTools)
library(Matrix)
library(ggplot2)
library(fields)
library(stringr)
library(seriation)
library(AUCell)
library(RColorBrewer)
library(pheatmap)

setwd('/projects/sle_jul_23_gabibov/')


source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

subsets_cols <- subsets_cols_level2

date <- '14Jul'

data_folder <-paste0('/projects/sle_jul_23_gabibov/luad/gene_modules/data/', date, '/')

path_to_intgr_seurat <- "/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds"

integrated <- readRDS(path_to_intgr_seurat)

integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel2

data_matrix <- readRDS(paste0("/projects/sle_jul_23_gabibov/luad/gene_modules/data/", date, "/data_matrix.rds"))
cormatrix <- readRDS(paste0("/projects/sle_jul_23_gabibov/luad/gene_modules/data/", date, "/cor_matrix.rds"))
modules <- readRDS(paste0("/projects/sle_jul_23_gabibov/luad/gene_modules/data/", date, "/modules.rds"))

# Filter out modules with less than 3 genes
#modules <- Filter(function(x) length(x) >= 3, modules)

genes <- unlist(modules, use.names = F)

nmods <- length(modules)

#создаем матрицу модули х гены, все по нулям
mod_mat <- matrix(0, nrow = nmods, ncol = length(genes), dimnames = list(1:nmods, genes))

#для каждого модуля гены, кот содерж именно в нем, становятся 1, т е теперь марица из нулей и 0, каждый ген встречается лишь в одном модуле
for(iter in 1:length(modules)){
  mod_mat[iter, modules[[iter]]] <- 1
}

cells <- colnames(data_matrix)

#получаем матрицу модули х клетки, с учетом UMI (очень грубое объяснение - для каждого модуля показано какой экспрессионный объем занимают гены из него для каждой клетки)
#mod_scores <- t(t(mod_mat%*%crc_matrix)/colSums(crc_matrix)) кажется что не надо делать (в ориг скрипте делали) т к у меня уже нормализ каунты в crc_matrix
mod_scores <- mod_mat%*%data_matrix[genes,]

mod_scores_df <- as.data.frame(t(mod_scores))

colnames(mod_scores_df) <- paste0("AvgModule", 1:50)

mod_scores_df$barcode_join <- rownames(mod_scores_df)

#делаем порядок генов в таблице матрицы корреляций генов таким чтобы гены из одного модуля были рядом
cmat <- cormatrix[genes, genes]

avg_cor <- array(NA, nmods, dimnames = list(names(modules)))

#в avg_cor записываю средние корреляции между генами раздельно по модулям
for(mod_iter in names(avg_cor)){
  genes <- modules[[mod_iter]]
  mat.iter <- cmat[genes, genes]
  diag(mat.iter) <- NA
  avg_cor[mod_iter] <- mean(mat.iter, na.rm = T)
}

#в списке с модулями и генами гены упорядочиваю гены так, что гены, коррелирующие с остальными генами из модуля, находятся первее
for(mod_iter in names(modules)){
  genes <- modules[[mod_iter]]
  mat.iter <- cmat[genes, genes]
  diag(mat.iter) <- NA
  modules[[mod_iter]] <- genes[order(rowSums(mat.iter, na.rm = T), decreasing = T)]
}

modules5top <- lapply(modules, function(x) head(x, 5))

# Add scores to object
integrated <- AddModuleScore(
  integrated,
  features = modules,
  name = "GeneModule"
)


# Add scores to object
# integrated <- AddModuleScore(
#   integrated,
#   features = modules5top,
#   name = "Gene5topModule"
# )

# # Run AUCell
# 
# integrated_kras <- subset(integrated, subset = KRAS != 'Unknown')
# 
# names(modules) <- paste0("AUCell", seq_along(modules))
# 
# exprMatrix <- GetAssayData(integrated_kras, assay = "RNA", slot = "data")
# 
# geneSets <- lapply(modules, unique)
# geneSets <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
# 
# cellsAUC <- AUCell_calcAUC(modules, geneSets)
# 
# aucMatrix <- as.data.frame(getAUC(cellsAUC))
# 
# integrated_kras <- AddMetaData(integrated_kras, t(aucMatrix))


DotPlot(integrated, features = integrated@meta.data %>% select(contains('GeneModule')) %>% names(), group.by = 'CelltypeAnnotationLow') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(data_folder, "gene_modules.pdf"), bg="white", height = 8, width = 40, dpi = 300, limitsize = FALSE)


DotPlot(integrated, features = integrated@meta.data %>% select(contains('Gene5topModule')) %>% names(), group.by = 'CelltypeAnnotationLow') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(data_folder, "gene5top_modules.pdf"), bg="white", height = 8, width = 40, dpi = 300, limitsize = FALSE)

DotPlot(integrated, features = integrated@meta.data %>% select(contains('GeneModule')) %>% names(), group.by = 'seurat_clusters') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(data_folder, "gene_modules_seuratClusters.pdf"), bg="white", height = 8, width = 40, dpi = 300, limitsize = FALSE)


write.csv(integrated@meta.data, paste0(data_folder, 'meta_50modules.csv'))


modules_df <- do.call(rbind, lapply(names(modules), function(name) {
  data.frame(module = name, gene = modules[[name]], stringsAsFactors = FALSE)
}))

write.csv(modules_df, paste0(data_folder, 'modules_df.csv'), row.names = FALSE)

#---------------------------------------------------------------FIND INTERESTING MODULES

modules_scores <- read.csv(paste0(data_folder, 'meta_50modules.csv'))

modules_scores <- modules_scores %>% filter(KRAS != 'Unknown')

relevant_modules <- modules_scores %>% select(contains('GeneModule')) %>% names()

relevant_modules

results_list <- list()

interest_celltypes <- c('Tem Th1 like', 
                        'Temra',
                        'Tem Th1 like GZMK',
                        'Tem Th1 like Exhausted PD1 high',
                        'Tem Th17 like CCR6 high',
                        'Tem Th17 like CCR6 low ',
                        'Tem Tfh like CXCL13 high',
                        'Tem Tfh like CXCL13 low',
                        'T reg CD25 high',
                        'T reg CD25 low',
                        'Exhausted PD1 high')


for (ct in unique(modules_scores$CelltypeAnnotationLow)) {
  
  meta_subset <- modules_scores[modules_scores$CelltypeAnnotationLow == ct, ]
  
  res_ct <- data.frame(
    CellType = ct,
    Module = relevant_modules,
    Median_WT = NA,
    Median_mut = NA,
    Diff = NA,
    P_value = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(relevant_modules)) {
    module <- relevant_modules[i]
    
    scores_wt <- meta_subset[meta_subset$KRAS == "WT", module]
    scores_mut <- meta_subset[meta_subset$KRAS == "mut", module]
    
    if (length(scores_wt) == 0 || length(scores_mut) == 0) {
      next
    }
    
    test <- wilcox.test(scores_wt, scores_mut)
    
    res_ct$Median_WT[i]  <- median(scores_wt, na.rm = TRUE)
    res_ct$Median_mut[i] <- median(scores_mut, na.rm = TRUE)
    res_ct$Diff[i]       <- res_ct$Median_mut[i] - res_ct$Median_WT[i]
    res_ct$P_value[i]    <- test$p.value
  }
  
  results_list[[ct]] <- res_ct
}

final_results <- do.call(rbind, results_list)

final_results_sorted <- final_results[order(-abs(final_results$Diff)), ]

final_results_sorted <- na.omit(final_results_sorted)

write.csv(final_results_sorted, paste0(data_folder, 'modules_stats.csv'), row.names = FALSE)



Celltype <- 'Tem Tfh like CXCL13 high'
Celltype <- 'Temra'
Celltype <- 'T reg CD25 high'
Celltype <- 'Tem Th1 like Exhausted PD1 high'
Celltype <- 'Trm Th1 like'
Celltype <- 'Cycling high'
Celltype <- 'IFN-induced'
Celltype <- 'TcmII'
Celltype <- 'Tem Th17 like CCR6 low'

ggplot(modules_scores %>% filter(CelltypeAnnotationLow == Celltype), aes(x = KRAS, y = GeneModule11, fill = KRAS)) +
  geom_boxplot()

#---------------------------------------------------------------SAVE GOOD MODULES

modules_scores <- read.csv(paste0(data_folder, 'meta_50modules.csv'))

bad_modules <- c(6, 12, 14, 17, 20, 21, 40, 42, 43, 45, 49, 50)

bad_modules <- paste0('GeneModule', bad_modules)

integrated <- AddMetaData(integrated, modules_scores %>% tibble::column_to_rownames(var = 'X'))

integrated@meta.data$CelltypeAnnotationLow <- factor(integrated@meta.data$CelltypeAnnotationLow, levels = c('TcmI',
                                                                                                            'TcmII',
                                                                                                            'Trm',
                                                                                                            'Tcm Th1-Th17 like',
                                                                                                            'Trm Th17 like CCR6 low',
                                                                                                            'Tem Th17 like CCR6 high',
                                                                                                            'Tem Th17 like CCR6 low',
                                                                                                            'Trm Th1 like',
                                                                                                            'Tem Th1 like',
                                                                                                            'Tem Th1 like GZMK',
                                                                                                            'Tem Th1 like Exhausted PD1 high',
                                                                                                            'Temra',
                                                                                                            'Tcm Tfh like CXCL13 low',
                                                                                                            'Trm Tfh like CXCL13 low',
                                                                                                            'Tem Tfh like CXCL13 high',
                                                                                                            'Tem Tfh like CXCL13 low',
                                                                                                            'Exhausted PD1 high',
                                                                                                            'T reg CD25 high',
                                                                                                            'T reg CD25 low',
                                                                                                            'IFN-induced',
                                                                                                            'Cycling high',
                                                                                                            'Cycling low',
                                                                                                            'Artefact',
                                                                                                            'Unknown'))

integrated <- subset(integrated, subset = CelltypeAnnotationLow %in% setdiff(unique(integrated@meta.data$CelltypeAnnotationLow), c('Artefact', 'Unknown')))

good_modules <- setdiff(integrated@meta.data %>% select(contains('GeneModule')) %>% names(), bad_modules)

DotPlot(integrated, features = good_modules, group.by = 'subclusters_final') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(data_folder, "gene_modules_subclusters.png"), bg="white", height = 20, width = 20, dpi = 300, limitsize = FALSE)

Celltype <- 'Tem Tfh like CXCL13 high'
Celltype <- 'Tem Th1 like Exhausted PD1 high'
Celltype <- 'T reg CD25 high'
Celltype <- 'Temra'

# 4E GENE MODULE 2------------------------------------------------------------

Celltypes <- c('Tem Tfh like CXCR5 high CXCL13 high', 'PD1 high', 'Tem Tfh like CXCR5 low CXCL13 high', 'Cytotoxic PD1 high')

modules_scores <- modules_scores %>% filter(KRAS != 'Unknown')

# Define colors
color_mut <- col_vector[19]
color_wt <- col_vector[30]

# Compute medians per group and cell type
medians <- modules_scores %>%
  filter(CelltypeAnnotationLow %in% Celltypes) %>%
  group_by(CelltypeAnnotationLow, KRAS) %>%
  summarise(med = median(GeneModule2, na.rm = TRUE), .groups = "drop")

# Prepare data for lines
medians_lines <- medians %>%
  mutate(color = ifelse(KRAS == "KRASmut", color_mut, color_wt))

# Compute median differences for center labels
median_diff <- medians %>%
  pivot_wider(names_from = KRAS, values_from = med) %>%
  mutate(diff = mut - WT)

# Plot
ggplot(modules_scores %>% filter(CelltypeAnnotationLow %in% Celltypes), 
       aes(x = KRAS, y = GeneModule2, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  facet_wrap(~CelltypeAnnotationLow, ncol = 2) +
  scale_fill_manual(values = c("mut" = color_mut, "WT" = color_wt)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2") +
  
  # Add dashed median lines
  geom_hline(data = medians_lines,
             aes(yintercept = med, color = KRAS),
             linetype = "dashed", size = 0.6, show.legend = FALSE) +
  scale_color_manual(values = c("mut" = color_mut, "WT" = color_wt)) +
  
  # Add text showing difference of medians
  geom_text(data = median_diff,
            aes(x = 1.5, y = max(modules_scores$GeneModule2, na.rm = TRUE) * 0.65,
                label = paste0("Δ = ", round(diff, 2))),
            inherit.aes = FALSE, size = 4, color = "black") +
  
  # Add p-values
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    label.x = 1.5,
    size = 5,
    color = 'red'
  )

ggsave(paste0(data_folder, "4E_geneMod2.pdf"), bg="white", height = 5, width = 7, dpi = 300, limitsize = FALSE)

#---------------------------------------------------------------GENE MODULE 2 PATIENT LEVEL

module2_patient <- modules_scores %>%
  group_by(Patient, CelltypeAnnotationLow, KRAS) %>%
  summarise(GeneModule2 = median(GeneModule2), .groups = "drop")

# Step 1: Filter to relevant cell types
filtered_cells <- modules_scores %>%
  filter(CelltypeAnnotationLow %in% Celltypes)

# Step 2: Count cells per (Patient, CelltypeAnnotationLow)
cell_counts <- filtered_cells %>%
  group_by(Patient, CelltypeAnnotationLow) %>%
  summarise(n_cells = n(), .groups = "drop")

# Step 3: Keep only (Patient, Celltype) pairs with at least 20 cells
valid_pairs <- cell_counts %>%
  filter(n_cells >= 20)

# Step 4: Filter the summarized module2_patient accordingly
module2_patient_filtered <- module2_patient %>%
  semi_join(valid_pairs, by = c("Patient", "CelltypeAnnotationLow")) %>%
  filter(CelltypeAnnotationLow %in% Celltypes)

module2_patient_filtered <- module2_patient_filtered %>% filter(GeneModule2 >= 0)



ggplot(module2_patient_filtered %>% filter(CelltypeAnnotationLow %in% Celltypes), 
       aes(x = KRAS, y = GeneModule2, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  facet_wrap(~CelltypeAnnotationLow, ncol = 2) +
  scale_fill_manual(values = c("mut" = color_mut, "WT" = color_wt)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2") +
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    label.x = 1.5,
    size = 5,
    color = 'red'
  )

#---------------------------------------------------------------GENE MODULE 1------------------------------------------------------------

Celltypes <- c('Tem Th1 like GZMK high', 'Tem Th1 like GZMA high', 'Cytotoxic PD1 high', 'Cytotoxic TEMRA')

modules_scores <- modules_scores %>% filter(KRAS != 'Unknown')

ggplot(modules_scores %>% filter(CelltypeAnnotationLow %in% Celltypes), aes(x = KRAS, y = GeneModule5, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  facet_wrap('CelltypeAnnotationLow', ncol = 2) +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2") +
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    label.x = 1.5,
    size = 5,
    color = 'red'
  )







modules_scores_pretreatment <- modules_scores %>% filter(is.na(Treatment))

ggplot(modules_scores_pretreatment %>% filter(CelltypeAnnotationLow %in% Celltypes), aes(x = KRAS, y = GeneModule2, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  facet_wrap('CelltypeAnnotationLow', ncol = 2, scales = 'free') +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2") +
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    label.x = 1.5,
    size = 5,
    color = 'red'
  )

ggsave(paste0(data_folder, "geneMod2_preT.png"), bg="white", height = 5, width = 7, dpi = 300, limitsize = FALSE)

#---------------------------------------------------------------EXHAUSTION GENE MODULE------------------------------------------------------------

exhaustion_geneset <- list(c('PDCD1', 'TOX', 'TOX2', 'TIGIT', 'CTLA4', 'HAVCR2', 'IKZF2', 'ENTPD1', 'NR4A1', 'NR4A2', 'NR4A3'))

integrated <- AddModuleScore(
  integrated,
  features = exhaustion_geneset,
  name = "ExhaustionGeneset"
)

Celltypes <- c('Tem Tfh like CXCL13 high', 'Tem Th1 like Exhausted PD1 high', 'Exhausted PD1 high')

metadata <- integrated@meta.data %>% filter(KRAS != 'Unknown')

ggplot(metadata %>% filter(CelltypeAnnotationLow %in% Celltypes), aes(x = KRAS, y = ExhaustionGeneset1, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  facet_wrap('CelltypeAnnotationLow', ncol = 2) +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2") +
  ggpubr::stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    label.x = 1.5,
    size = 5,
    color = 'red'
  )

#---------------------------------------------------------------HEATMAP------------------------------------------------------------

modules_scores <- read.csv(paste0(data_folder, 'meta_50modules.csv'))

modules_scores <- filter(modules_scores, CelltypeAnnotationLow %in% setdiff(unique(modules_scores$CelltypeAnnotationLow), c('Unknown', 'Artefact')))

modules_scores$CelltypeAnnotationLow <- ifelse(modules_scores$CelltypeAnnotationLow == 'TcmI', 'Tcm/Naive', modules_scores$CelltypeAnnotationLow)
modules_scores$CelltypeAnnotationLow <- ifelse(modules_scores$CelltypeAnnotationLow == 'TcmII', 'Tcm/Naive', modules_scores$CelltypeAnnotationLow)

celltypes_order <- c('Tcm/Naive',
                     'Trm',
                     'Tcm Th1-Th17 like',
                     'Trm Th17 like CCR6 low',
                     'Tem Th17 like CCR6 high',
                     'Tem Th17 like CCR6 low',
                     'Trm Th1 like',
                     'Tem Th1 like',
                     'Tem Th1 like GZMK',
                     'Tem Th1 like Exhausted PD1 high',
                     'Temra',
                     'Tcm Tfh like CXCL13 low',
                     'Trm Tfh like CXCL13 low',
                     'Tem Tfh like CXCL13 high',
                     'Tem Tfh like CXCL13 low',
                     'Exhausted PD1 high',
                     'T reg CD25 high',
                     'T reg CD25 low',
                     'IFN-induced',
                     'Cycling high',
                     'Cycling low')

bad_modules <- c(6, 12, 14, 17, 20, 21, 40, 42, 43, 45, 49, 50, 28, 18, 4, 29, 13, 33, 48, 16, 32, 44, 24, 47, 19, 25)

bad_modules <- paste0('GeneModule', bad_modules)

good_modules <- setdiff(colnames(modules_scores %>% select(contains('GeneModule'))), bad_modules)

good_modules_scores <- modules_scores[, good_modules]

good_modules_scores$Subset <- modules_scores$CelltypeAnnotationLow

good_modules_scores_mean <- good_modules_scores %>%
  group_by(Subset) %>%
  summarise(across(starts_with("GeneModule"), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames('Subset')

good_modules_scores_mean <- good_modules_scores_mean[celltypes_order, ]

good_modules_scores_scaled <- scale(good_modules_scores_mean)

p1 <- pheatmap(good_modules_scores_scaled,
               cluster_rows = FALSE,
               cluster_cols = TRUE,
               #annotation_row = annotation_row,
               #annotation_col = annotation_col,
               #annotation_colors = c(annotation_colors_col), # Add if needed - annotation_row_colors
               fontsize_row = 8,
               fontsize_col = 8,
               #annotation_legend = TRUE,
               main = "Heatmap of Gene Modules by Cell Type",
               color = colorRampPalette(c("blue", "white", "#FDDA0D", "red"))(100),
               breaks = seq(-2, 4, length.out = 101),
               border_color = 'white')

p1


png(file=paste0(data_folder, "heatmap_genemodules.png"), width=20, height=15, units="cm", res=300)
draw(p1)
dev.off()


Celltype <- 'Tem Tfh like CXCL13 high'
Celltype <- 'Tem Th1 like Exhausted PD1 high'

modules_scores <- modules_scores %>% filter(KRAS != 'Unknown')

ggplot(modules_scores %>% filter(CelltypeAnnotationLow == Celltype), aes(x = KRAS, y = GeneModule2, fill = KRAS)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  geom_violin(trim = FALSE, width = 0.5, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Module Score: GeneModule2")



#---------------------------S9?------------------------------------MODULE FROM THE ARTICLE------------------------------------------------------------

neoantigen_genes <- c('MT-ND1', 'MALAT1', 'VIM', 'EMP3', 'CXCL13', 'NR3C1', 'ADGRG1', 'NMB', 'ITM2A', 'ETV7', 'COTL1', 'B2M', 'IGFL2')

neoantigen_modules <- list(c('CXCL13', 'NR3C1', 'ADGRG1', 'NMB', 'ITM2A', 'ETV7', 'COTL1', 'B2M', 'IGFL2'))

integrated <- AddModuleScore(
  integrated,
  features = neoantigen_modules,
  name = "NeoModule"
)

integrated_clean <- subset(integrated, subset = CelltypeAnnotationLow %in% setdiff(unique(integrated@meta.data$CelltypeAnnotationLow), c('Artefact', 'Unknown')))

# Get average expression matrix
avg_exp <- AverageExpression(integrated_clean, features = neoantigen_genes, group.by = "CelltypeAnnotationLow", slot = "data")$RNA

celltype_order <- intersect(names(subsets_cols), colnames(avg_exp))
gene_order <- intersect(neoantigen_genes, rownames(avg_exp))

# Subset and reorder matrix
ordered_mat <- avg_exp[gene_order, celltype_order]

# Scale per gene (optional, improves contrast)
ordered_mat_scaled <- scale(t(ordered_mat))

# ANNOTATION CELLS
annotation_row <- data.frame(CellType = factor(celltype_order, levels = celltype_order))
rownames(annotation_row) <- celltype_order

annotation_row_colors <- list(CellType = subsets_cols)

# ANNOTATION GENES

annotation_row <- data.frame(Regulation = c(rep("Down-regulated", 4), rep("Up-regulated", 9)))
rownames(annotation_row) <- neoantigen_genes

# Annotation color
ann_colors <- list(
  Regulation = c("Down-regulated" = "#4075D7", "Up-regulated" = "#CF1C1C")
)


p1 <- pheatmap(as.matrix(t(ordered_mat_scaled)),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_row = annotation_row,
               #annotation_col = annotation_col,
               annotation_colors = c(ann_colors), # Add if needed - annotation_row_colors
               fontsize_row = 10,
               fontsize_col = 12,
               annotation_legend = TRUE,
               main = "Heatmap of Marker Genes by Cell Type",
               color = colorRampPalette(c("blue", "white", "#FDDA0D", "red"))(100),
               breaks = seq(-2, 4, length.out = 101),
               border_color = 'white')


png(file=paste0(data_folder, "heatmap_neoantigen_genes.png"), width=18,height=14,units="cm", res=600)
grid::grid.newpage()
grid::grid.draw(p1$gtable)
#draw(p1)
dev.off()


ggplot(integrated_clean@meta.data %>% filter(KRAS != 'Unknown'), aes(x = CelltypeAnnotationLow, y = NeoModule1, fill = KRAS)) +
  #geom_point(position = position_jitterdodge(), size = 0.5, alpha = 0.4) +
  #geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Enrichment score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        panel.grid.major.x = element_blank()) +
  geom_hline(yintercept = 0.75, color = 'red', linetype = 'dashed')

ggsave(paste0(data_folder, "neoantigen_module.png"), bg="white", height = 5, width = 10, dpi = 300, limitsize = FALSE)

# 2B EXHAUSTION AND CYTOTOXICITY MODULES------------------------------------------------------------

celltypes <- c('PD1 high', 'Cytotoxic PD1 high', 'Cytotoxic TEMRA') #'Tem Tfh like', 'Tem Th1 like', 

subsets_cols <- subsets_cols_level1

integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel1

exhaustion_module <- list(c('PDCD1', 'HAVCR2', 'CTLA4', 'LAG3', 'TIGIT'))

cytotoxicity_module <- list(c('GZMB', 'GZMA', 'GZMH', 'GZMK', 'PRF1', 'GNLY', 'NKG7'))

integrated <- AddModuleScore(
  integrated,
  features = exhaustion_module,
  name = "ExhaustionModule"
)

integrated <- AddModuleScore(
  integrated,
  features = cytotoxicity_module,
  name = "CytotoxicityModule"
)

integrated@meta.data %>% select(c('CelltypeAnnotationLevel1', 'CelltypeAnnotationLevel2', 'Patient', 'KRAS', 'orig.ident', 'ExhaustionModule1', 'CytotoxicityModule1')) %>%
  write.csv(paste0(data_folder, 'cytotox_and_exhaust_modules.csv'), row.names = TRUE)



meta_cyt_exh <- read.csv(paste0(data_folder, 'cytotox_and_exhaust_modules.csv'), row.names = 1)

meta_cyt_exh <- meta_cyt_exh %>% filter(CelltypeAnnotationLevel1 %in% celltypes)

meta_cyt_exh$ExhaustionModule1_norm <- (meta_cyt_exh$ExhaustionModule1 - min(meta_cyt_exh$ExhaustionModule1)) / (max(meta_cyt_exh$ExhaustionModule1) - min(meta_cyt_exh$ExhaustionModule1))
meta_cyt_exh$CytotoxicityModule1_norm <- (meta_cyt_exh$CytotoxicityModule1 - min(meta_cyt_exh$CytotoxicityModule1)) / (max(meta_cyt_exh$CytotoxicityModule1) - min(meta_cyt_exh$CytotoxicityModule1))


meta_cyt_exh$ExhaustionModule1_norm <- meta_cyt_exh$ExhaustionModule1
meta_cyt_exh$CytotoxicityModule1_norm <- meta_cyt_exh$CytotoxicityModule1



summary_df <- meta_cyt_exh %>%
  group_by(CelltypeAnnotationLevel1) %>%
  summarise(
    median_exh = median(ExhaustionModule1_norm, na.rm = TRUE),
    q1_exh = quantile(ExhaustionModule1_norm, 0.25, na.rm = TRUE),
    q3_exh = quantile(ExhaustionModule1_norm, 0.75, na.rm = TRUE),
    median_cyt = median(CytotoxicityModule1_norm, na.rm = TRUE),
    q1_cyt = quantile(CytotoxicityModule1_norm, 0.25, na.rm = TRUE),
    q3_cyt = quantile(CytotoxicityModule1_norm, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- summary_df %>%
  tidyr::pivot_longer(
    cols = -CelltypeAnnotationLevel1,
    names_to = c("stat", "score"),
    names_sep = "_"
  ) %>%
  tidyr::pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  mutate(score = dplyr::recode(score,
                               exh = "Inhibitory Receptor Score",
                               cyt = "Cytotoxicity Score"))

plot_df$CelltypeAnnotationLevel1 <- factor(plot_df$CelltypeAnnotationLevel1, levels = c('PD1 high', 'Cytotoxic PD1 high', 'Tem Th1 like', 'Cytotoxic TEMRA'))

ggplot(plot_df, aes(x = CelltypeAnnotationLevel1, y = median)) +
  geom_point(aes(color = score), size = 3) +
  geom_errorbar(aes(ymin = q1, ymax = q3, color = score), width = 0.1) +
  geom_line(aes(group = 1, color = score)) +
  facet_wrap(~score, scales = "free_y", ncol = 1) +  # stacked vertically
  scale_color_manual(values = c("Inhibitory Receptor Score" = "#339B0E", "Cytotoxicity Score" = "#D42B13")) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = "none"
  ) #+
#ylim(0, 0.8)
#scale_y_continuous(limits = c(0, 0.8), breaks = c(0, 0.8))

ggsave(paste0(data_folder, "2B_cyt_and_exh_scores.pdf"), bg="white", height = 4, width = 4, dpi = 600, limitsize = FALSE)


# S8 Clonality green and blue ------------------------------------------------


library(Seurat)
library(tidyverse)
library(here)
library(pals)
library(ggpubr)
library(RColorBrewer)

setwd('/projects/sle_jul_23_gabibov/')


source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

subsets_cols <- subsets_cols_level2

# Annotations should be stored under metadata Subset column, TCR cdr3 beta under TCR2nt,
# donor under patient_id

figures_folder <- "/home/dlukyanov/All_R_projects/projects/scRNA_seq/our/CD4_projection/figures/nslc_kras"
data_folder <- '/projects/sle_jul_23_gabibov/luad/tcr_analysis/data/'

path_to_intgr_seurat <- "/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds"

intgr_full <- read_rds(path_to_intgr_seurat)


integrated_backup <- intgr_full

intgr_full@meta.data$CelltypeAnnotationLow <- intgr_full@meta.data$CelltypeAnnotationLevel2

#-----------------------------------------------------------------DELETE UNKNOWN AND ARTEFACT CELL TYPES

intgr_full <- subset(intgr_full, subset = CelltypeAnnotationLevel2  != 'Unknown')

#intgr <- intgr_full
intgr <- subset(intgr_full, subset = KRAS != 'Unknown')


intgr@meta.data$patient_id <- intgr@meta.data$Patient
intgr@meta.data$Subset <- intgr@meta.data$CelltypeAnnotationLow
intgr@meta.data$TCR2nt <- intgr@meta.data$cdr3_nt2

intgr$Barcode <- row.names(intgr[[]])

#  Delete cells without TCR data
cells_keep <- rownames(intgr@meta.data)[!is.na(intgr@meta.data$TCR2nt)]
intgr <- subset(intgr, cells = cells_keep)

#  No less than N cells in each Subset in each patient

N <- 50
intgr_filtered_cells <- group_by(intgr[[]], Subset, patient_id) %>% filter(n() > N)
intgr <- subset(intgr, cells = intgr_filtered_cells$Barcode)

intgr$Barcode <- NULL
intgr$barcode <- NULL # in case there is such column

# extract UMAP coordinates for cells 

umap_tx <- intgr@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(intgr[[]]) %>%
  rownames_to_column(var = "barcode")

# clonality plot  

# 100% - all TRB clones

# umap_tx %>%
#   dplyr::select(patient_id, Subset, TCR2nt) %>%
#   drop_na() %>% # we take only TRB containing cells %>%
#   group_by(patient_id, Subset) %>%
#   mutate(n_cells_in_cluster = n()) %>%
#   distinct() %>% # keep only unique clonotypes
#   group_by(patient_id, Subset,  n_cells_in_cluster) %>%
#   summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
#   group_by(patient_id) %>%
#   mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
#          fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
#          n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

umap_tx %>%
  dplyr::select(patient_id, Subset, TCR2nt, KRAS) %>%
  drop_na() %>% # we take only TRB containing cells %>%
  group_by(patient_id, Subset, KRAS) %>%
  mutate(n_cells_in_cluster = n()) %>%
  distinct() %>% # keep only unique clonotypes
  group_by(patient_id, Subset, KRAS, n_cells_in_cluster) %>%
  summarize(n_clones_in_cluster = n(), .groups = "drop") %>%
  group_by(patient_id, KRAS) %>%
  mutate(n_total_clones_in_donor = sum(n_clones_in_cluster),
         fraction_of_TRB_clones_in_cluster_within_all_TRB_clones = n_clones_in_cluster / n_total_clones_in_donor,
         n_total_cells_in_donor = sum(n_cells_in_cluster)) -> clonality_data

clonality_data %>%
  mutate(ratio = n_clones_in_cluster / n_cells_in_cluster) %>%
  mutate(Subset = factor(Subset, ordered = T,
                         levels=names(subsets_cols))) -> clonality_data

clonality_data$colors <- subsets_cols[clonality_data$Subset]

#write.csv(clonality_data, paste0(data_folder, 'clonality_df_10cellsF.csv'))

# Filter out cell types with less than 2 patients UNIQUE FEATURE

K <- 0

patients_per_celltypes <- clonality_data %>%
  distinct(patient_id, KRAS, Subset) %>%
  group_by(Subset, KRAS) %>%
  summarise(n_patients = n(), .groups = "drop") %>%
  ungroup %>%
  group_by(Subset) %>%
  filter(n_distinct(KRAS) == 2) %>%
  filter(all(n_patients[KRAS == "WT"] >= K),
         all(n_patients[KRAS == "mut"] >= K)) %>%
  ungroup()

clonality_data <- clonality_data %>%
  filter(Subset %in% unique(patients_per_celltypes$Subset))

# Filter out cell types with less than 2 patients UNIQUE FEATURE

# FILTER BY MEDIAN ABSOLUTE DEVIATION

stats <- clonality_data %>%
  group_by(Subset, KRAS) %>%
  summarise(
    med = median(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones),
    mad_val = mad(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones),
    .groups = "drop"
  )

clonality_data <- clonality_data %>%
  left_join(stats, by = c("Subset", "KRAS")) %>%
  mutate(
    mad_score = abs(fraction_of_TRB_clones_in_cluster_within_all_TRB_clones - med) / mad_val
  )

# PLOT DISTRIBUTION OF MAD SCORES

threshold <- quantile(clonality_data$mad_score, 0.9, na.rm = TRUE)

mad_score_distribution <- ggplot(clonality_data, aes(x = mad_score)) +
  geom_histogram(bins = 105) +
  theme_bw() +
  xlab('MAD score') +
  ylab('Count') +
  geom_vline(xintercept = threshold, color = 'red', linetype = 'dashed')

mad_score_distribution

#ggsave(plot = mad_score_distribution, paste0(figures_folder, "madScore_prepostT_09perc_madF_min50cells.png"), bg="white", height = 8, width = 10, dpi = 300, units = "cm")


# FILTER OUT ONLY SAMPLES WITH HIGH MAD

clonality_data <- clonality_data %>%
  filter(mad_score <= threshold)

#write.csv(clonality_data, paste0(data_folder, 'clonality_df_10cellsF_madF.csv'))

# FILTER BY MEDIAN ABSOLUTE DEVIATION

clonality_data$clonality <- 1 - clonality_data$ratio 

plotClonalityGroups <- function(x){
  # Shows clonality between groups inside covariate
  clonality_data <- x
  c1 <- ggplot(clonality_data, aes(x = Subset, y = clonality, fill = KRAS))  +
    coord_flip() +
    geom_boxplot(alpha = 0.7) +
    # stat_compare_means(method = "anova", show.legend = FALSE, 
    #                    label.x = length(unique(clonality_data$Subset)), 
    #                    label.y = 0.5,
    #                    size = 3) +
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = paste0(max(clonality_data$Subset)),
    #                    label.y = 0, hide.ns = T) +
    geom_point(key_glyph = "point",position=position_dodge(width = .5)) +
    theme_minimal() +
    ggtitle("TCR\u03b2 clonality") +
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1))) +
    scale_color_identity(guide = "legend", labels = names(subsets_cols), breaks = subsets_cols) +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          panel.grid.major.y = element_blank()) +
    #scale_x_log10() +
    #scale_y_log10() +
    xlab("") +
    ylab("Clonality")  + 
    NoLegend() + 
    #ylim(0, 1) +
    scale_fill_manual(values = col_vector[c(19,30)]) +
    scale_y_reverse(breaks = seq(0, 1, 0.25))
  
  c2 <- ggplot(data = clonality_data, aes(x = Subset, fill = KRAS)) +
    coord_flip() +
    geom_bar(stat = 'count', position = 'dodge') +
    xlab('') +
    ylab('') +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    #scale_y_continuous(
    #  limits = c(0, y_limit),
    #  breaks = seq(0, y_limit, by = 5)) +
    scale_fill_manual(values = col_vector[c(19,30)])
  
  fig <- ggpubr::ggarrange(c1, c2, nrow = 1, align = 'h', common.legend = TRUE, legend = 'bottom', widths = c(5, 1), vjust = 0.1)
  
  plot(fig)
}

plot_clonality <- plotClonalityGroups(clonality_data)

ggsave(plot = plot_clonality, paste0(figures_folder, "/S8_clonality_kras_prepostT_NomadF_min50cells.pdf"), bg="white", height = 7, width = 10, dpi = 600)


# S9A S9B Neoantigen module -------------------------------------------------------

library(Seurat)
library(dplyr)
library(DescTools)
library(Matrix)
library(ggplot2)
library(fields)
library(stringr)
library(seriation)
library(AUCell)
library(RColorBrewer)
library(pheatmap)

setwd('/projects/sle_jul_23_gabibov/')


source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

subsets_cols <- subsets_cols_level1

date <- '14Jul'

figures_folder <-paste0('/projects/sle_jul_23_gabibov/luad/gene_modules/neoantigen_figs/')

path_to_intgr_seurat <- "/projects/sle_jul_23_gabibov/luad/annotation_big/rds/luad_rpca_annotated_21Jun25.rds"

integrated <- readRDS(path_to_intgr_seurat)

integrated@meta.data$CelltypeAnnotationLow <- integrated@meta.data$CelltypeAnnotationLevel1

#---------------------------------------------------------------MODULE FROM THE ARTICLE

neoantigen_genes <- c('MT-ND1', 'MALAT1', 'VIM', 'EMP3', 'CXCL13', 'NR3C1', 'ADGRG1', 'NMB', 'ITM2A', 'ETV7', 'COTL1', 'B2M', 'IGFL2')

neoantigen_modules <- list(c('CXCL13', 'NR3C1', 'ADGRG1', 'NMB', 'ITM2A', 'ETV7', 'COTL1', 'B2M', 'IGFL2'))

integrated <- AddModuleScore(
  integrated,
  features = neoantigen_modules,
  name = "NeoModule"
)

integrated_clean <- subset(integrated, subset = CelltypeAnnotationLow != 'Unknown')

# Get average expression matrix
avg_exp <- AverageExpression(integrated_clean, assays = 'RNA', features = neoantigen_genes, group.by = "CelltypeAnnotationLow", slot = "data")$RNA

celltype_order <- intersect(names(subsets_cols), colnames(avg_exp))
gene_order <- intersect(neoantigen_genes, rownames(avg_exp))

# Subset and reorder matrix
ordered_mat <- avg_exp[gene_order, celltype_order]

# Scale per gene (optional, improves contrast)
ordered_mat_scaled <- scale(t(ordered_mat))

# ANNOTATION CELLS
annotation_row <- data.frame(CellType = factor(celltype_order, levels = celltype_order))
rownames(annotation_row) <- celltype_order

annotation_row_colors <- list(CellType = subsets_cols)

# ANNOTATION GENES

annotation_row <- data.frame(Regulation = c(rep("Down-regulated", 4), rep("Up-regulated", 9)))
rownames(annotation_row) <- neoantigen_genes

# Annotation color
ann_colors <- list(
  Regulation = c("Down-regulated" = "#4075D7", "Up-regulated" = "#CF1C1C")
)


p1 <- pheatmap(as.matrix(t(ordered_mat_scaled)),
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_row = annotation_row,
               #annotation_col = annotation_col,
               annotation_colors = c(ann_colors), # Add if needed - annotation_row_colors
               fontsize_row = 10,
               fontsize_col = 12,
               annotation_legend = TRUE,
               main = "Heatmap of Marker Genes by Cell Type",
               color = colorRampPalette(c("blue", "white", "#FDDA0D", "red"))(100),
               breaks = seq(-2, 4, length.out = 101),
               border_color = 'white')


pdf(file=paste0(figures_folder, "/S9A_heatmap_neoantigen_genes.pdf"), width=18 / 2.54,height=14 / 2.54)
draw(p1)
dev.off()


ggplot(integrated_clean@meta.data %>% filter(KRAS != 'Unknown'), aes(x = CelltypeAnnotationLow, y = NeoModule1, fill = KRAS)) +
  #geom_point(position = position_jitterdodge(), size = 0.5, alpha = 0.4) +
  #geom_violin(trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  scale_fill_manual(values = col_vector[c(19,30)]) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    panel.grid.major.x = element_blank()
  ) +
  ylab("Enrichment score") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        panel.grid.major.x = element_blank()) +
  geom_hline(yintercept = 0.75, color = 'red', linetype = 'dashed')

ggsave(paste0(figures_folder, "/S9B_neoantigen_module.pdf"), bg="white", height = 5, width = 10, dpi = 600, limitsize = FALSE)










