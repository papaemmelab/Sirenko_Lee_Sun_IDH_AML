library(dplyr)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
"%ni%" = Negate("%in%")
options(ggrepel.max.overlaps = Inf)
output_dir <- "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/manuscript/figure5/"

all.combined.new_filtered <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/IDH_AML_integrated_seurat_obejct.rds")

################################################################################################################################################
#### DEG wrapper function 
################################################################################################################################################

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




celltypes = c("MPP", "GMP", "CD14+ Mono")

## AML-08
donor.id = "AML-08"
seurat_object <- all.combined.new_filtered[, all.combined.new_filtered$donor.id %in% c(donor.id, paste0("Control-", stringr::str_pad(c(5:12), width = 1 ,pad = "0")))]
seurat_object <- seurat_object[, seurat_object$manual_annotation_v3_l3_final %in% celltypes]
seurat_object$stage_v2 <- seurat_object$stage_detail
seurat_object$stage_v2[grep(seurat_object$stage_v2, pattern = "CR")] <- "CR"

comparison = list(c("Dx", "healthy"), c("R1", "healthy"), c("R2", "healthy"), c("P1", "healthy"), list(c("R1", "R2", "P1"), "Dx"))
deg_list_AML08 <- DEGTest(seurat_object, 
                          stage_col= "stage_v2", 
                          celltype_col= "manual_annotation_v3_l3_final", 
                          comparison = comparison, 
                          celltypes= celltypes, 
                          logfc.threshold =0.1,
                          min.pct = 0.1)
## AML-07
donor.id = "AML-07"
seurat_object <- all.combined.new_filtered[, all.combined.new_filtered$donor.id %in% c(donor.id, paste0("Control-", stringr::str_pad(c(5:12), width = 1 ,pad = "0")))]
seurat_object <- seurat_object[, seurat_object$manual_annotation_v3_l3_final %in% celltypes]
comparison = list(c("Dx", "healthy"), c("R1", "healthy"), c("R2", "healthy"), list(c("R1", "R2"), "Dx"))
deg_list_AML07 <- DEGTest(seurat_object, 
                          stage_col= "stage_detail", 
                          celltype_col= "manual_annotation_v3_l3_final", 
                          comparison = comparison, 
                          celltypes= celltypes, 
                          logfc.threshold =0.1,
                          min.pct = 0.1)
## AML-10
donor.id = "AML-10"
seurat_object <- all.combined.new_filtered[, all.combined.new_filtered$donor.id %in% c(donor.id, paste0("Control-", stringr::str_pad(c(5:12), width = 1 ,pad = "0")))]
seurat_object <- seurat_object[, seurat_object$manual_annotation_v3_l3_final %in% celltypes]
comparison = list(c("Dx", "healthy"), c("R", "healthy"), c("R", "Dx"))
deg_list_AML10 <- DEGTest(seurat_object, 
                          stage_col= "stage_detail", 
                          celltype_col= "manual_annotation_v3_l3_final", 
                          comparison = comparison, 
                          celltypes= celltypes, 
                          logfc.threshold =0.1,
                          min.pct = 0.1)



#### filter the genes and reformat the data frames 
deg_list_filtered <- deg_list_AML10
  
for(i in 1:length(deg_list_filtered)){
  for(j in 1:length(deg_list_filtered[[1]])){
    deg_list_filtered[[i]][[j]]$gene <- rownames(deg_list_filtered[[i]][[j]])
    deg_list_filtered[[i]][[j]]$comparison <- names(deg_list_filtered[[i]][j])
    deg_list_filtered[[i]][[j]]$celltype <- names(deg_list_filtered[i])
    deg_list_filtered[[i]][[j]] <- subset(deg_list_filtered[[i]][[j]], p_val_adj < 0.01)
    
    gene_filtered <- deg_list_filtered[[i]][[j]]$gene[-grep(deg_list_filtered[[i]][[j]]$gene, pattern = "^MT-|^RPS|^RPL|XIST|^HBB|^HBA|^HBD|^HBM")]
    if (length(gene_filtered)>0){
      deg_list_filtered[[i]][[j]] <- subset(deg_list_filtered[[i]][[j]], gene %in% gene_filtered)
    }
  }
  # deg_write <- Reduce(rbind, deg_list_filtered[[i]])
  # write.table(deg_write, file = paste0(output_dir, names( deg_list_filtered[i]) ,".tsv"), quote = F, sep = "\t", row.names = F)
}

## select the top genes for the heatmap 

n_deg =20
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

  


### heatmap 
celltype_col= "manual_annotation_v3_l3_final"

DefaultAssay(seurat_object) <- "RNA"
celltypes = c("MPP", "GMP", "CD14+ Mono")
# celltypes="GMP"
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
  data_df<- cbind(data_df, celltype = seurat_object@meta.data[cells,c("stage_simple")])
  mean_df <- data_df %>% group_by(celltype) %>% summarize_all(mean)
  mean_df_scaled <- scale(mean_df[,-1])
  rownames(mean_df_scaled) <- mean_df$celltype
  attr(mean_df_scaled, "scaled:center") <- NULL
  attr(mean_df_scaled, "scaled:scale") <- NULL
  
  pdf(paste0(output_dir, donor.id, "_", celltypes[i], "-top_" ,n_deg, "_heatmap2_confirmed.pdf"), width = 2.5, height = sqrt(length(gene_interest))-0.5)
  p1 <- Heatmap(t(mean_df_scaled)[gene_interest,], name="Z-score", 
                #row_order = gene_interest,
                row_names_gp = gpar(fontsize = 8), 
                column_names_gp = gpar(fontsize = 8), 
                column_order = rownames(mean_df_scaled))
  print(p1)
  dev.off()
}




  
  
