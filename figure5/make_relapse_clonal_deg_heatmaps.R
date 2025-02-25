
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(ggpattern)
library(ggalluvial)
library(cowplot) 

"%ni%" = Negate("%in%")

source("/gpfs/home/sirenm01/sirenm01/projects/ref/plotting_helpers.R")

# load in latest seurat object 
int_obj_all = readRDS( "objects/idh_aml_object_all.rds")

pw = 5
ph = 5

# DEGs across multiple conditions
DEGTest <- function(seurat_object, stage_col, celltype_col, celltypes,  comparison, logfc.threshold =0.1, min.pct = 0.1){
  deg_list <- list()
  comparison_string <- unlist(lapply(comparison, function(x) paste0(paste0(x[[1]], collapse = ""), sep = "_", paste0(x[[2]], collapse = ""))))
  DefaultAssay(seurat_object) <- "RNA"
  for (i in 1:length(celltypes)){
    deg_single_celltype <- list()
    for (j in 1:length(comparison)){
      message("Run FindMarkers: ", comparison_string[j], " for ", celltypes[i])
      cells.compared1 <- colnames(seurat_object)[seurat_object[[stage_col]][[1]] == comparison[[j]][[1]] & seurat_object[[celltype_col]][[1]] %in% celltypes[i]]
      cells.compared2 <- colnames(seurat_object)[seurat_object[[stage_col]][[1]] == comparison[[j]][[2]] & seurat_object[[celltype_col]][[1]] %in% celltypes[i]]
      if(length(cells.compared1)>0 & length(cells.compared2)>0){
        deg <- FindMarkers(seurat_object, ident.1 = cells.compared1, ident.2 = cells.compared2, logfc.threshold = logfc.threshold, min.pct = min.pct)
        deg_single_celltype[[j]] <- deg
      } else {
        warning(paste0("No cells are present for the comparison: ", comparison_string[j], ". Skip the comparison"))
      }
    }
    names(deg_single_celltype) <- comparison_string
    deg_list[[i]] <- deg_single_celltype
  }
  names(deg_list) <- celltypes
  return(deg_list) 
}


########################
####### AML-17 #########
########################

## subset seurat object for cells from AML-17
meta_blasts = int_obj_all[[]] %>% dplyr::filter(manual_annotation_v3_l1_final %in% c("CD14+ Mono",    "MPP","cDC","GMP", "Progenitor_DC" ,   "Ery", "MEP","HSC","CD16+ Mono")) %>% 
dplyr::filter(donor.id == "AML-17")
meta_blasts$celltype = meta_blasts$manual_annotation_v3_l1_final

## assign cells to clonal/disease stage groups for DGE
meta_df = meta_blasts %>% 
  mutate(del7_group = case_when(
    (numbat == "chr7_del" & orig.ident %in% c("AML-17-D1", "AML-17-D2")) ~ "Dx-del7",
    (numbat == "chr7_del" & orig.ident %in% c("AML-17-R1", "AML-17-R2")) ~ "Rel-del7",
    (numbat == "normal" & orig.ident %in% c("AML-17-D1", "AML-17-D2")) ~ "Dx-nk",
    (numbat == "normal" & orig.ident %in% c("AML-17-R1", "AML-17-R2")) ~ "Rel-nk",
    TRUE ~ "Other"
  )) %>% dplyr::select(donor.id, del7_group)

AML_17 = subset(int_obj_all, cells = rownames(meta_df))
AML_17 = AddMetaData(AML_17, meta_df)

# Focus on cells belonging to clonal/disease stage groups 
AML_17 = subset(AML_17, del7_group %in% c("Rel-del7", "Dx-del7","Rel-nk", "Dx-nk"))

celltypes = c("MPP") ## focus on MPP

# run DGE
comparison = list(c("Rel-del7", "Dx-del7"), c("Rel-nk", "Dx-nk"), c("Rel-del7","Rel-nk"), c("Dx-del7","Dx-nk"))
deg_list_AML17 <- DEGTest(AML_17, 
                          stage_col= "del7_group", 
                          celltype_col= "predicted_idh_celltype", 
                          comparison = comparison, 
                          celltypes= celltypes, 
                          logfc.threshold =0,
                          min.pct = 0.1)
# saveRDS(deg_list_AML17, "degs/deg_list_AML17_fc0_del7_clones_addDx_Dx.rds")

deg_list_AML17 = readRDS("degs/deg_list_AML17_fc0_del7_clones_addDx_Dx.rds")

deg_list_filtered <- deg_list_AML17
seurat_object = AML_17
 donor.id = "AML17_del7_clones"
  
for(i in 1:length(deg_list_filtered)){
  for(j in 1:length(deg_list_filtered[[1]])){
    deg_list_filtered[[i]][[j]]$gene <- rownames(deg_list_filtered[[i]][[j]])
    deg_list_filtered[[i]][[j]]$comparison <- names(deg_list_filtered[[i]][j])
    deg_list_filtered[[i]][[j]]$celltype <- names(deg_list_filtered[i])
    deg_list_filtered[[i]][[j]] <- subset(deg_list_filtered[[i]][[j]], p_val_adj < 0.01)
    
    gene_filtered <- deg_list_filtered[[i]][[j]]$gene[-grep(deg_list_filtered[[i]][[j]]$gene, pattern = "^MT-|^RPS|^RPL|XIST|^HBB|^HBA|^HBD|^HBM|^MALAT")]
    if (length(gene_filtered)>0){
      deg_list_filtered[[i]][[j]] <- subset(deg_list_filtered[[i]][[j]], gene %in% gene_filtered)
    }
  }
  # deg_write <- Reduce(rbind, deg_list_filtered[[i]])
  # write.table(deg_write, file = paste0(output_dir, names( deg_list_filtered[i]) ,".tsv"), quote = F, sep = "\t", row.names = F)
}


## select the top genes for the heatmap 
n_deg =15
p_val_cutoff = 1e-2
deg_list_selected <- list()

for (i in 1:length(deg_list_filtered)){
  deg_list_selected_temp <- list()
  for(j in 1:length(deg_list_filtered[[1]])){
    deg_selected <- deg_list_filtered[[i]][[j]] %>% mutate(DEG = ifelse(avg_log2FC>0, "up", "down"))
    deg_selected <- deg_selected %>% subset(p_val_adj < p_val_cutoff) %>% group_by(DEG) %>% top_n(n_deg, wt = abs(avg_log2FC)) %>% arrange(desc(avg_log2FC))
    deg_list_selected_temp[[j]] <- deg_selected
  }
  deg_list_selected[[i]] <- deg_list_selected_temp
  names(deg_list_selected[[i]]) <-  names(deg_list_filtered[[i]])
}
names(deg_list_selected) <- names(deg_list_filtered)


## annotate chr7 genes 
library(AnnotationDbi)
library(org.Hs.eg.db)

# Mark genes located on chromosome 7 with an asterisk
get_chr7_genes <- function(genes) {
  gene_info <- mapIds(org.Hs.eg.db, keys = genes, column = "CHR", keytype = "SYMBOL")
  chr7_genes <- genes[gene_info == "7"]
  # genes[genes %in% chr7_genes] <- paste0(genes[genes %in% chr7_genes], "*")
  return(chr7_genes)
}


### make heatmap 
celltype_col= "predicted_idh_celltype"
DefaultAssay(seurat_object) <- "RNA"
celltypes = c("MPP")

for (i in 1:length(celltypes)){
    cells <- colnames(seurat_object)[seurat_object[[celltype_col]][[1]] %in% celltypes[i]]
    # gene_interest <- unique(c(deg_list_selected[[celltypes[i]]]$deg_relapse$gene, 
    #                           deg_list_selected[[celltypes[i]]]$deg_diagnosis$gene,
    #                           deg_list_selected[[celltypes[i]]]$deg_relapse_diagnosis$gene))
    # 
    gene_interest <- unique(unlist(lapply(deg_list_selected$MPP, function(x) x$gene)))

    ## remove some noncoding RNAs
    if(length(grep(gene_interest, pattern = "^AC[0-9][0-9][0-9][0-9][0-9]"))>0){
        gene_interest <- gene_interest[-grep(gene_interest, pattern = "^AC[0-9][0-9][0-9][0-9][0-9]")]
    }
    data_df <- FetchData(seurat_object, vars = gene_interest, slot = "data", cells = cells)
    data_df<- cbind(data_df, del7_group = seurat_object@meta.data[cells,c("del7_group")])
    data_df$del7_group = factor(data_df$del7_group, levels = c("Dx-nk", "Dx-del7", "Rel-nk", "Rel-del7"))
    mean_df <- data_df %>% group_by(del7_group) %>% summarize_all(mean)
    mean_df_scaled <- scale(mean_df[,-1])
    rownames(mean_df_scaled) <- mean_df$del7_group
    attr(mean_df_scaled, "scaled:center") <- NULL
    attr(mean_df_scaled, "scaled:scale") <- NULL
  

    # Remove columns with NA/NaN/Inf values
    mean_df_scaled <- mean_df_scaled[, !colSums(is.na(mean_df_scaled) | is.nan(mean_df_scaled) | is.infinite(mean_df_scaled))]

    # Apply the chromosome 7 marker
    chromosome_7_genes <- get_chr7_genes(gene_interest)
    colnames(mean_df_scaled) <- ifelse(colnames(mean_df_scaled) %in% chromosome_7_genes, paste0("*", colnames(mean_df_scaled) ), colnames(mean_df_scaled))

    gene_interest0 = colnames(mean_df_scaled)

    #####   horizontal heatmaps w/ genes on the bottom 
    pdf(paste0(output_dir, donor.id, "_", celltypes[i], "-top_" ,n_deg, "_heatmap_del7_chr7_genes_addDx_Dx_20250108_horizonal.pdf"), height = 1.75,  width= sqrt(length(gene_interest)))
    p1 <- Heatmap((mean_df_scaled), name="Z-score", 
                    row_order = c("Dx-nk"  ,  "Dx-del7" , "Rel-nk" ,  "Rel-del7"),
                    show_row_dend = FALSE,
                    row_names_gp = gpar(fontsize = 6), 
                    column_names_gp = gpar(fontsize = 8) 
                    # column_order = colnames(mean_df_scaled)
                    )
    print(p1)
    dev.off()

}