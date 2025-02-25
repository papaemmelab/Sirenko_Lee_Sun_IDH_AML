
### instead make a scillus heatmap with the DEGs from IDH and TET2 CH + normals 

## findallmarkers 
setwd("/gpfs/data/aifantislab/home/sirenm01/projects/idh/aml16_17_18_19_20_21_22_CCUS/")
source("/gpfs/home/sirenm01/sirenm01/projects/ref/plotting_helpers.R")

ch_mye = readRDS("objects/ch_mye_tet2_idh_ccus.rds")
ch_all = readRDS( "objects/ch_all_tet2_idh_ccus.rds")
options(max.overlaps = 80)
options(ggrepel.max.overlaps = 40)


setwd("/gpfs/data/aifantislab/home/sirenm01/projects/idh/revisions_relapse/")
source("/gpfs/home/sirenm01/sirenm01/projects/ref/plotting_helpers.R")
output_dir = "/gpfs/data/aifantislab/home/sirenm01/projects/idh/revisions_relapse/"




## DEGS
idh_tet2_pb_df = ch_mye[[]] %>% dplyr::filter(orig.ident %in% c("CH-01-T1-PB-P4","CH-01-T1-PB-P7","CH-05-T1-PB-LD",  "N-01-T1-PB-L",  "N-01-T1-PB-LD", "N-01-T1-PB-Stringent","N-02-T1-PB-L", "TET2_CCUS_08","TET2_CCUS_10", "TET2-05", "TET2-06", "TET2-07"), cell_type %in% c("CD14+ Mono"))
idh_tet2_pb = subset(ch_mye, cells = rownames(idh_tet2_pb_df)) 

idh2_ch = idh_tet2_pb_df %>% dplyr::filter(orig.ident %in% c("CH-01-T1-PB-P4","CH-01-T1-PB-P7","CH-05-T1-PB-LD"))
tet2_ch = idh_tet2_pb_df %>% dplyr::filter(orig.ident %in% c( "TET2_CCUS_08","TET2_CCUS_10",  "TET2-05", "TET2-06", "TET2-07"))
normal_pb = idh_tet2_pb_df %>% dplyr::filter(orig.ident %in% c( "N-01-T1-PB-L",  "N-01-T1-PB-LD", "N-01-T1-PB-Stringent","N-02-T1-PB-L"))
Idents(idh_tet2_pb) = "orig.ident"
# c6_obj_small <- subset(c6_obj, downsample = 70) 

# Add a metadata column based on genotype group
idh_tet2_pb$group <- case_when(
  idh_tet2_pb$orig.ident %in% c("CH-01-T1-PB-P4", "CH-01-T1-PB-P7", "CH-05-T1-PB-LD") ~ "idh2_ch",
  idh_tet2_pb$orig.ident %in% c("TET2_CCUS_08", "TET2_CCUS_10", "TET2-05", "TET2-06", "TET2-07") ~ "tet2_ch",
  idh_tet2_pb$orig.ident %in% c("N-01-T1-PB-L", "N-01-T1-PB-LD", "N-01-T1-PB-Stringent", "N-02-T1-PB-L") ~ "normal_pb",
  TRUE ~ NA_character_ 
)

# SDownsample to 1000 cells per donor
# Use donor.id for idh2_ch and normal_pb; use orig.ident for tet2_ch
downsampled_cells <- idh_tet2_pb@meta.data %>% tibble::rownames_to_column('cb') %>%
  mutate(donor_id = case_when(
    group == "idh2_ch" | group == "normal_pb" ~ donor.id,
    group == "tet2_ch" ~ orig.ident,
    TRUE ~ NA_character_
  )) %>%
  group_by( donor_id) %>%
  sample_n(size = min(1000, n()), replace = FALSE) %>%
  dplyr::select('cb')

idh_tet2_pb_small <- subset(idh_tet2_pb, cells = downsampled_cells$cb)

idh_tet2_pb_small@meta.data = idh_tet2_pb_small@meta.data %>% 
  mutate(donor_id = case_when(
    group == "idh2_ch" | group == "normal_pb" ~ donor.id,
    group == "tet2_ch" ~ orig.ident,
    TRUE ~ NA_character_
  )) 

# Set group as the identity class
Idents(idh_tet2_pb_small) <- "group"

DefaultAssay(idh_tet2_pb_small) = "RNA"


# Find DEGs with FindAllMarkers
markers <- FindAllMarkers(
  object = idh_tet2_pb_small,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0
)


label_genes = unique(c(
  "NFKBIA", "TNFAIP3", "CXCL8", "HIF1A", "ICAM1", "SERPINA1", # Inflammation
  "SPI1", "CEBPB", "CEBPD", "IRF1",                         # Myeloid Differentiation
  "LYZ", "PTPRC", "VCAN",                                   # Monocyte Function
  "HLA-DQB1", "ICAM1", "TMEM176B", "TMEM176A",              # Antigen Presentation
  "JUN", "HSP90AA1", "MAPK14", "NR4A1", "FOSB",              # Proliferation
  "FTH1", "EIF1", "CD99", "TAGLN2", "IFITM3", "SERPINA1", "DDT", "HSPA5", 
  "PIM3", "HIF1A", "NFKBIA", "NR4A1", "NR4A2",
  "SPI1", "GADD45B", "LYZ", 
  "ICAM1", "IRF1", "MAPK14", "PTPRC", 
  "EGR1", "CEBPD", "IER2", "FOSB", "CXCR4", "FOS", "IER3", 
  "IFNGR2", "EIF1", "NINJ1", "BCL3", "G0S2", "KLF10", "BCL2A1", "SMAD3", 
  "NFKBIE", "BCL6", "NFKB2", "REL", "ICAM1", "HIF1A", "PIM1", "DDIT4", "HMOX1", 
  "ANXA2", "BHLHE40", "HSPA5", "ZFP36", "CDKN1A", "TIMP1", "IL10RA", "IRF1", "TNFRSF1B", 
  "IFITM1", "IFITM2", "IFITM3", "CD44", "MAP3K8", "TGFB1", "ITGA4", "CD99", "SAT1"
))

# Select top DEGs to plot
top_markers0 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30) %>%
  ungroup()

already_included = label_genes[(label_genes %in% top_markers0$gene)]

top_selected <- markers %>% 
  dplyr::filter((gene %in% label_genes)) %>% 
  dplyr::filter(!(gene %in% already_included)) 

top_markers = rbind(top_markers0, top_selected) 

# Display the top markers
print(top_markers, n =30)

# Define the donor IDs
donor_ids <- c(
  "CH-01", "CH-02", 
  "healthy_MSK_01", "healthy_MSK_02", 
  "TET2-01", "TET2-02", "TET2-03", 
  "TET2-04", "TET2-05"
)

# Assign colors to donor IDs
donor_colors <- c(
  # Dark purples for CH-01 and CH-05
  "CH-01" = "#6b4294",  # Indigo
  "CH-02" = "#9071af",  # Royal Purple

  # Beiges for healthy_MSK
  "healthy_MSK_01" = "#e3ba8a" , ##F5DEB3",  # Wheat
  "healthy_MSK_02" = "#f1ddc4", ##FFE4C4",  # Bisque

  # Shades of red/orange for TET2 samples
  "TET2-01" = "#b2d8d8",  # Orange Red
  "TET2-02" = "#66b2b2",  # Tomato
  "TET2-03" = "#008080",  # Coral
  "TET2-04" = "#006666",  # Light Salmon
  "TET2-05" = "#004c4c"   # Salmon
)

library(Scillus)
library(ComplexHeatmap)
library(circlize)

# saveRDS(markers, "degs_findallmarkers_CH_idh_tet2_pb_logFC0_downsamp1000.rds")
# write.csv(markers, "degs_findallmarkers_CH_idh_tet2_pb_logFC0_downsamp1000.csv")

markers = readRDS("degs_findallmarkers_CH_idh_tet2_pb_logFC0_downsamp1000.rds")


idh_tet2_pb_small = ScaleData(idh_tet2_pb_small, features = rownames(idh_tet2_pb_small))


idh_tet2_pb_small$group <- factor(
  idh_tet2_pb_small$group,
  levels = c("normal_pb", "idh2_ch", "tet2_ch")
)




plot_heatmap_with_labels_right <- function(dataset, 
                                           markers,
                                           sort_var = c('seurat_clusters'),
                                           n = 8, 
                                           anno_var, 
                                           anno_colors,
                                           anno_label, 
                                           hm_limit = c(-2, 0, 2), 
                                           hm_colors = c("#4575b4","white","#d73027"),
                                           row_font_size = 12,
                                           labeled_genes) {
        
        # Extract scaled data
        mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
        
        # Filter and arrange genes
        if (is.data.frame(markers)) {
            genes <- get_top_genes(dataset, markers, n)
        } else if (is.character(markers)) {
            genes <- markers
        } else {
            stop('Incorrect input of markers')
        }
        
        mat <- mat[match(genes, rownames(mat)),]
        
        # Create annotation metadata
        anno <- dataset@meta.data %>%
                tibble::rownames_to_column(var = "barcode") %>%
                arrange(!!!syms(sort_var))
        
        # Reorder the matrix by annotations
        mat <- t(mat)
        mat <- mat[match(anno$barcode, rownames(mat)),]
        mat <- t(mat)
    
        # Generate annotations
        annos <- list()
        for (i in seq_along(1:length(anno_var))) {
            value <- anno[[anno_var[i]]]
            
            if (is.numeric(value)) {
                col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
                                      anno_colors[[i]])
                ha <- HeatmapAnnotation(a = value,
                                        col = list(a = col_fun),
                                        border = TRUE,
                                        annotation_label = anno_label[i], 
                                        annotation_legend_param = list(title = anno_label[i])
                                        )
            } else {
                col <- list(a = anno_colors[[i]]) #setNames(anno_colors[[i]], levels(factor(value))))
                ha <- HeatmapAnnotation(a = value,
                                        col = col,
                                        border = TRUE,
                                        annotation_label = anno_label[i],
                                        annotation_legend_param = list(title = anno_label[i])
                                        )
            }
            annos[[i]] <- ha
        }
        
        annos <- do.call(c, annos)
        annos@gap <- rep(unit(1, "mm"), length(annos))
        
        # # Create row annotations for gene labels
        # gene_labels <- rowAnnotation(
        #     genes = 
        # Create row annotations for gene labels
        gene_labels <- rowAnnotation(
            genes = anno_mark(
                at = which(rownames(mat) %in% labeled_genes),
                labels = rownames(mat)[rownames(mat) %in% labeled_genes], 
                labels_gp = gpar(fontsize = row_font_size), 
                link_width = unit(9, "mm")
            )
        )
        # (
        #         at = which(rownames(mat) %in% labeled_genes),
        #         labels = rownames(mat)[rownames(mat) %in% labeled_genes]
        #     )
        # )
        
        # Build heatmap
        ht <- Heatmap(mat,
                      cluster_rows = TRUE,
                      cluster_columns = FALSE,
                      heatmap_legend_param = list(direction = "vertical",
                                                  legend_width = unit(6, "cm"),
                                                  title = "Expression"),
                      col = colorRamp2(hm_limit, hm_colors),
                      show_column_names = FALSE,
                      row_names_side = "left", 
                      row_names_gp = gpar(fontsize = row_font_size),
                      right_annotation = gene_labels, # Add gene annotations
                      top_annotation = annos,
                      column_split = anno[[anno_var[1]]],
                #       row_split = rep(c("A", "B"), each = length(markers) / length(unique(anno[[anno_var[1]]]))), 
                      gap = unit(1, "mm")
                      )
        
        draw(ht, 
             heatmap_legend_side = "right",
             annotation_legend_side = "right")
}

set_colors <- function(pal, n) {
        
        if (all(pal %in% rownames(brewer.pal.info))) {
                num <- c()
                for (i in seq(length(pal))) {
                        num[i] <- brewer.pal.info[pal[i],][[1]]
                }
                full_pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
        } else if (all(are_colors(pal))) {
                full_pal <- pal
        } else {
                stop('Incorrect palette setup. Please input valid RColorBrewer palette names or color names.')
        }
                
        if (n <= length(full_pal)) {
                return(full_pal[1:n])
        } else {
                warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.", 
                        immediate. = TRUE)
                return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
        }
}




p = plot_heatmap_with_labels_right(idh_tet2_pb_small, 
                         markers = unique(top_markers$gene),
                         labeled_genes = label_genes, 
                         sort_var = c('group', 'donor_id'),
                        #  n = 8, 
                         anno_var = c('group', 'donor_id'), 
                         anno_label = c("Genotype", "Donor ID"),
                         anno_colors = list(
    group = c("idh2_ch" = "#2b2134", "tet2_ch" = "#003535", "normal_pb" = "#f6e9d9"), 
    donor_id = donor_colors
  ), 
                         hm_limit = c(-2, 0, 2), 
                         hm_colors = c(scales::muted("#4575b4"),"white",scales::muted("#d73027")),
                         row_font_size = 8) 
pdf("degs/de_figures/heatmap_subset_idh_tet2_pb_small_labeled_130.pdf", width = 8, height = 10)
print(p)
dev.off()
dev.off()








