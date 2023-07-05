library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(ggrepel)
output_dir <- "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/manuscript/figure1/"


all.combined.new_filtered <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_adt_processed_101322.rds")

DefaultAssay(all.combined.new_filtered) <- "RNA"
all.combined.new_filtered$manual_annotation_v3_l3_final <- factor(all.combined.new_filtered$manual_annotation_v3_l3_final, 
                                                                  levels =  c("HSC", "MPP", "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage",
                                                                              "Progenitor_DC", "cDC", "MEP", "Ery", "Platelet", "Plasma", "Progenitor_B", "B",
                                                                              "pDC", "CD4+ T", "CD8+ T", "NK", "Stroma"))

Idents(all.combined.new_filtered) <- all.combined.new_filtered$manual_annotation_v3_l3_final

####################################################################################################
##### RNA MARKERS ##################################################################################
####################################################################################################
rna_markers <- FindAllMarkers(all.combined.new_filtered, only.pos = T, logfc.threshold = 1, min.pct = 0.2, max.cells.per.ident = 5e+3)


### RNA - top 30 genes from each cell type 
gene_interests <- rna_markers %>% subset(p_val_adj<1e-200) %>% group_by(cluster) %>% top_n(wt= avg_log2FC, n=30) %>% data.frame() 
gene_interests <- gene_interests %>% group_by(gene) %>% filter(avg_log2FC == max(avg_log2FC))
gene_interests <- unique(gene_interests$gene)
data_df <- FetchData(all.combined.new_filtered, vars = gene_interests, slot = "data")
data_df<- cbind(data_df, 
                stage = all.combined.new_filtered$stage_simple, 
                celltype = all.combined.new_filtered@meta.data[,c("manual_annotation_v3_l3_final")])

library(ComplexHeatmap)
gene_labels = c("AVP", "SPINK2", "CRHBP","SOCS2", "SMIM24", "C1QTNF4", "CDK6", "SNHG29", "CD34", "NPM1", "TRH", "AZU1", "MPO", "ELANE", "S100A9",
                "S100A8", "VCAN", "LYZ", "CD14", "CEBPD", "FCGR3A", "CDKN1C",  "FCER1G", "C1QC", "C1QA", "C1QB", "TUBA1B", "FCER1A",
                "HLA-DQA1" ,  "HLA-DPA1",   "HLA-DPB1",   "HLA-DRA", "CD1C", "TPSB2", "GATA2", "SOX4" ,"CD82", "HBB","HBA2", "HBA1",
                "PPBP", "PF4", "TUBB1", "IGKC", "IGHA1", "IGHG1", "IGHM", "IGHD", "PAX5","CD79A", "MS4A1", "CD22", "IRF7", "IRF8", "IL7R", "LTB", "BCL11B", "IL32", 
                "GZMH", "GZMB", "NKG7", "KLRB1",  "CXCL12", "APOE", "IGFBP5", "GAS6", "TF", "VCAM1" ,"COL1A2")


### all cells 
mean_df <- data_df %>% group_by(celltype) %>% 
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  dplyr::select(where(is.numeric)) %>% as.matrix()
rownames(mean_df) <- levels(data_df$celltype)
mean_df_scaled <- t(scale(mean_df))

pdf(paste0(output_dir, "All_combined_RNA_heatmap.pdf"), width = 5, height = 12)
Heatmap(mean_df_scaled[gene_interests,], name="Z-score", 
        row_order = gene_interests,
        column_order = c("HSC", "MPP", "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage",
                         "Progenitor_DC", "cDC", "MEP", "Ery", "Platelet", "Plasma", "Progenitor_B", "B",
                         "pDC", "CD4+ T", "CD8+ T", "NK", "Stroma") ) + 
  rowAnnotation(link = anno_mark(at = match(gene_labels, gene_interests), 
                                 labels = gene_labels,
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

dev.off()


### only healthy cells 
mean_df <- data_df %>% 
  subset(stage == "healthy") %>% 
  group_by(celltype)  %>% 
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  dplyr::select(where(is.numeric)) %>% as.matrix()
rownames(mean_df) <- levels(data_df$celltype)
mean_df_scaled <- t(scale(mean_df))

pdf(paste0(output_dir, "healthy_cells_combined_RNA_heatmap.pdf"), width = 5, height = 12)
Heatmap(mean_df_scaled[gene_interests,], name="Z-score", 
        row_order = gene_interests,
        column_order = c("HSC", "MPP", "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage",
                         "Progenitor_DC", "cDC", "MEP", "Ery", "Platelet", "Plasma", "Progenitor_B", "B",
                         "pDC", "CD4+ T", "CD8+ T", "NK", "Stroma") ) + 
  rowAnnotation(link = anno_mark(at = match(gene_labels, gene_interests), 
                                 labels = gene_labels,
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

dev.off()

### without healthy cells 
mean_df <- data_df %>% 
  subset(stage != "healthy") %>% 
  group_by(celltype)  %>% 
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  dplyr::select(where(is.numeric)) %>% as.matrix()
rownames(mean_df) <- levels(data_df$celltype)
mean_df_scaled <- t(scale(mean_df))

pdf(paste0(output_dir, "without_healthy_cells_combined_RNA_heatmap.pdf"), width = 5, height = 12)
Heatmap(mean_df_scaled[gene_interests,], name="Z-score", 
        row_order = gene_interests,
        column_order = c("HSC", "MPP", "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage",
                         "Progenitor_DC", "cDC", "MEP", "Ery", "Platelet", "Plasma", "Progenitor_B", "B",
                         "pDC", "CD4+ T", "CD8+ T", "NK", "Stroma") ) + 
  rowAnnotation(link = anno_mark(at = match(gene_labels, gene_interests), 
                                 labels = gene_labels,
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

dev.off()


####################################################################################################
##### ADT MARKERS ##################################################################################
####################################################################################################
DefaultAssay(all.combined.new_filtered) <- "ADT"
## use datasets with the 270 antibody panel
all.combined.new_filtered_subset <-all.combined.new_filtered[, all.combined.new_filtered$batch_name %in% c("AML0024", "AML1371_1" , "AML1371_2", 
                                                                                                           "AML2123_1", "AML2123_2", "AML4340_1", "AML4340_2")]
all.combined.new_filtered_subset <- all.combined.new_filtered_subset[,colSums(all.combined.new_filtered_subset) >0]
all.combined.new_filtered_subset <- NormalizeData(all.combined.new_filtered_subset, normalization.method = 'CLR', margin = 2) %>% ScaleData()
adt_markers_clr <- FindAllMarkers(all.combined.new_filtered_subset, only.pos = T, logfc.threshold = 0.2, min.pct = 0.1)
adt_markers_clr %>% subset(p_val_adj<1e-20) %>% group_by(cluster) %>% top_n(wt= avg_log2FC, n=5) %>% data.frame() %>% dplyr::select(avg_log2FC, p_val_adj, cluster, gene)

gene_interests <- rownames(all.combined.new_filtered_subset)
data_df <- FetchData(all.combined.new_filtered_subset, vars = gene_interests, slot = "data")
data_df<- cbind(data_df, celltype = all.combined.new_filtered_subset@meta.data[,c("manual_annotation_v3_l3_final")])

mean_df <- data_df %>% group_by(celltype) %>% 
  summarise(across(where(is.numeric), ~ mean(.x))) %>% 
  dplyr::select(where(is.numeric)) %>% as.matrix()
rownames(mean_df) <- levels(data_df$celltype)
mean_df_scaled <- t(scale(mean_df))

## top 10 ADTs 
gene_interests <- adt_markers_clr %>% subset(p_val_adj<1e-20) %>% group_by(cluster) %>% top_n(wt= avg_log2FC, n=10) %>% data.frame() 
gene_interests <- gene_interests %>% group_by(gene) %>% filter(avg_log2FC == max(avg_log2FC))
gene_interests <- unique(gene_interests$gene)


pdf(paste0(output_dir, "All_combined_ADT_heatmap.pdf"), width = 5, height = 12)
Heatmap(mean_df_scaled[gene_interests,], name="Z-score", 
        row_order = gene_interests,
        column_order = c("HSC", "MPP", "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage",
                         "Progenitor_DC", "cDC", "MEP", "Ery", "Platelet", "Plasma", "Progenitor_B", "B",
                         "pDC", "CD4+ T", "CD8+ T", "NK", "Stroma") )
dev.off()




