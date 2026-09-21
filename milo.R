#-----------------------------------------------------------------LIBRARIES------------------------------------------------------
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


#-----------------------------------------------------------------HELPERS---------------------------------------------------------

source('luad/helpers/colors.R')
source('luad/helpers/functions.R')

#-----------------------------------------------------------------READ DATA------------------------------------------------------

print('Read data')

figures_folder <- '/path_to_folder/milo/figures/'
data_folder <- '/path_to_folder/milo/data/'

integrated <- readRDS('/path_to_folder/integrated.rds')

#-----------------------------------------------------------------DELETE UNKNOWN AND ARTEFACT CELL TYPES------------------------------------------------------

integrated <- subset(integrated, subset = CelltypeAnnotationLevel1 != 'Unknown')

#-----------------------------------------------------------------FILTER BY KRAS STATUS------------------------------------------------------

#integrated_kras <- subset(integrated, subset = KRAS != 'Unknown')

integrated_kras <- integrated

#-----------------------------------------------------------------CREATE MILO OBJECT------------------------------------------------------

print('Create Milo object')
#integrated[["RNA"]] <- as(integrated[["RNA"]], Class="Assay")

integrated_kras_sce <- as.SingleCellExperiment(integrated_kras, assay = 'RNA')

integrated_kras_milo <- Milo(integrated_kras_sce)

#-----------------------------------------------------------------CONSTRUCT KNN GRAPH------------------------------------------------------

print('Construct KNN graph')

integrated_kras_milo <- buildGraph(integrated_kras_milo, k = 30, d = 30, reduced.dim = 'INTEGRATED.RPCA')

#-----------------------------------------------------------------NHOODS------------------------------------------------------------------

print('Calculate Nhoods')

integrated_kras_milo <- makeNhoods(integrated_kras_milo, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = 'INTEGRATED.RPCA')

plotNhoodSizeHist(integrated_kras_milo)

cowplot::ggsave2('nhoods_size.png', path = figures_folder, height = 15, width = 20, units = 'cm', dpi = 300)

#-----------------------------------------------------------------COUNT CELLS IN NBHOODS------------------------------------------------------------

print('Count cells')

integrated_kras_milo <- countCells(integrated_kras_milo, meta.data = as.data.frame(colData(integrated_kras_milo)), sample="Patient")

#----------------------------------------------------------------DESIGN MATRIX------------------------------------------------------------

print('Create matrix')

design <- data.frame(colData(integrated_kras_milo))[,c("Patient", "KRAS", "orig.ident")]

## Convert batch info from integer to factor
design$orig.ident <- as.factor(design$orig.ident) 
design <- distinct(design)
rownames(design) <- design$Patient

design

#----------------------------------------------------------------NBHOOD CONNECTIVITY------------------------------------------------------------

print('Calculate connectivity')

integrated_kras_milo <- calcNhoodDistance(integrated_kras_milo, d=30, reduced.dim = "INTEGRATED.RPCA")

#----------------------------------------------------------------DIFFERENTIAL ABUNDANCE TESTING------------------------------------------------------------

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

#----------------------------------------------------------------SAVE MILO------------------------------------------------------------

print('Save results')

saveRDS(integrated_kras_milo, paste0(data_folder, 'integrated_milo_final.rds'))
write.csv(da_results, paste0(data_folder, 'da_milo_final.csv'))
write.csv(da_results_batch, paste0(data_folder, 'da_milo_final_batch.csv'))
write.csv(da_results_patient, paste0(data_folder, 'da_milo_final_patient.csv'))

#----------------------------------------------------------------READ MILO------------------------------------------------------------

# integrated_kras_milo <- readRDS(paste0(data_folder, 'integrated_milo_final.rds'))
# da_results <- read.csv(paste0(data_folder, 'da_milo_annotated.csv'))

#----------------------------------------------------------------ASSIGN NHOODS TO DISCRETE CLUSTERS------------------------------------------------------------

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)


da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel1")
da_results <- annotateNhoods(integrated_kras_milo, da_results, coldata_col = "CelltypeAnnotationLevel2")

da_results$CelltypeAnnotationLevel1 <- factor(da_results$CelltypeAnnotationLevel1, levels = names(subsets_cols_level1))
da_results$CelltypeAnnotationLevel2 <- factor(da_results$CelltypeAnnotationLevel2, levels = setdiff(names(subsets_cols_level2), 'Unknown'))

plotDAbeeswarm(da_results, group.by = "CelltypeAnnotationLevel2")

cowplot::ggsave2('da_plot_final.png', path = figures_folder, height = 20, width = 20, units = 'cm', dpi = 300)

#----------------------------------------------------------------FILTER NHOODS BY FREQUENCY IN CLUSTERS------------------------------------------------------------

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

cowplot::ggsave2('da_plot_final_filtered.png', path = figures_folder, height = 20, width = 25, units = 'cm', dpi = 300)

#----------------------------------------------------------------GOODNESS-OF-FIT TESTS------------------------------------------------------------

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

write.csv(metrics_summary, paste0(data_folder, 'metrics_summary_fdr005_freq09.csv'))

#----------------------------------------------------------------PLOT NHOODS ON KNN GRAPH------------------------------------------------------------

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

cowplot::ggsave2(filename = file.path(figures_folder, 'umaps_final_filtrered.png'),
                 plot = combined_plot2,
                 height = 30, width = 60, units = 'cm', dpi = 300)











#------------------------------------------------------------------------------------------------------------------------------------------------
