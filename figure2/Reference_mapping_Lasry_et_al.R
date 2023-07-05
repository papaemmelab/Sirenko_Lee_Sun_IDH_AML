library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(UCell)
library(ggpubr)
library(ComplexHeatmap)
library(UCell)
"%ni%" = Negate("%in%")
options(ggrepel.max.overlaps = Inf)

output_dir = "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/integrated/MSK_collaboration/diagnosis/"

### load the data and run SPCA

all.combined.new_filtered <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_081522.rds")
all.combined.new_filtered$manual_annotation_v3_l3_final <- as.character(all.combined.new_filtered$manual_annotation_v3_l1_final)
all.combined.new_filtered$manual_annotation_v3_l3_final[all.combined.new_filtered$manual_annotation_v3_l3_final=="Immature Mono"]  <- "MPP"

aml_dx.obj <- readRDS("/gpfs/data/aifantislab/home/nadorb01/share/seurat_obj_adult.rds")

DefaultAssay(all.combined.new_filtered) <- "RNA"
all.combined.new_filtered <- NormalizeData(all.combined.new_filtered)
all.combined.new_filtered <- FindVariableFeatures(all.combined.new_filtered, selection.method = "vst", nfeatures = 2000)
all.combined.new_filtered <- ScaleData(all.combined.new_filtered)
all.combined.new_filtered <- RunPCA(all.combined.new_filtered)
all.combined.new_filtered <- FindNeighbors(all.combined.new_filtered, graph.name= "RNA_snn")
all.combined.new_filtered <- RunSPCA(all.combined.new_filtered, assay = 'RNA', graph = 'RNA_snn')

all.combined.new_filtered <- FindNeighbors(
  object = all.combined.new_filtered,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

all.combined.new_filtered <- RunUMAP(all.combined.new_filtered, return.model =T, reduction = "integrated_pca", dims = 1:30, reduction.name = "rpca_umap", reduction.key = "RPCAUMAP_")

###################################################################################################################################
### find anchors between the reference dataset and the query dataset
DefaultAssay(aml_dx.obj) <- "RNA"

aml_dx.obj$donor.id <- aml_dx.obj$samples
aml_dx.obj$stage <- "diagnosis"
aml_dx.obj$stage[grep(aml_dx.obj$samples, pattern = "Control")] <- "healthy"
aml_dx.obj_subset <- aml_dx.obj[,aml_dx.obj$stage != "healthy"]
aml_dx.obj_subset <- aml_dx.obj_subset[, aml_dx.obj_subset$donor.id %ni% c(unique(all.combined.new_filtered$donor.id))]
aml_dx.obj.list <- SplitObject(aml_dx.obj_subset, split.by = "donor.id")

anchors <- list()
for (i in 1:length(aml_dx.obj.list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = all.combined.new_filtered,
    query = aml_dx.obj.list[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(aml_dx.obj.list)) {
  aml_dx.obj.list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = aml_dx.obj.list[[i]],
    reference = all.combined.new_filtered, 
    refdata = list(
      celltype = "manual_annotation_v3_l3_final"),
    reference.reduction = "spca",
    reduction.model = "rpca_umap"
  )
}

# Merge the batches 
aml_dx.obj_merged <- merge(aml_dx.obj.list[[1]], aml_dx.obj.list[2:length(aml_dx.obj.list)], merge.dr = c("ref.umap", "ref.spca"))
df <- aml_dx.obj_merged[['ref.umap']]@cell.embeddings
aml_dx.obj_merged[['rpca_umap']] <- CreateDimReducObject(embeddings= aml_dx.obj_merged[['ref.umap']]@cell.embeddings, key = "RPCAUMAP_")


### add gene module scores 
hallmark_gs <- read.delim("/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/Signatures/h.all.v7.5.1.symbols.gmt", header = F)
hallmark_genesets <- apply(hallmark_gs[,-c(1:2)], MARGIN = 1, function(x){unname(c(x))[x!=""]})
names(hallmark_genesets) <- hallmark_gs[,1]

HSC_Prog = c('SPINK2', 'ANGPT1', 'GUCY1A3', 'FAM30A', 'MMRN1', 'TPT1',
             'GAS5','RAB27B','TPM4','MSI2','GCSAML','SOCS2','EEF1A1', 'NRIP1','HOPX','CD34', 'TFPI',
             'TPSD1','PDZRN4','PCNP','PTPRCAP','FLT3', 'SMIM24', 'SELENOP','DAPK1','SMYD3',
             'ADGRG6', 'PIM1','MECOM', 'CEP70')

stemcell_genesets <- list(HSC_Prog)
names(stemcell_genesets) <- c("HSC_Prog")

genesets <- c(hallmark_genesets,stemcell_genesets)

all.combined.new_filtered@meta.data$IDH <- NULL
all.combined.new_filtered@meta.data$IDH[all.combined.new_filtered@meta.data$disease == "healthy"] <- "Control"
all.combined.new_filtered@meta.data$IDH[all.combined.new_filtered@meta.data$disease %in% c("AML", "CH")] <- "IDH"

col_remove = colnames(all.combined.new_filtered@meta.data)[grep(colnames(all.combined.new_filtered@meta.data), pattern = "_UCell")]
for (i in 1:length(col_remove)){
  all.combined.new_filtered[[col_remove[i]]] <- NULL
}

all.combined.new_filtered$id <- 'reference'
aml_dx.obj_merged$id <- 'query'
refquery <- merge(all.combined.new_filtered, aml_dx.obj_merged, merge.dr = c("rpca_umap"))
refquery[["spca"]] <- merge(all.combined.new_filtered[["spca"]], aml_dx.obj_merged[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50, 
                    reduction.name = "new_umap", reduction.key = "NEWUMAP_")


my.matrix <- GetAssayData(refquery, assay = "RNA", slot = "data")
scores <- ScoreSignatures_UCell(my.matrix, features=genesets)
refquery = AddMetaData(refquery, data.frame(scores))

### merge the celltype annotation column 
refquery@meta.data[colnames(aml_dx.obj_merged), "manual_annotation_v3_l3_final"] <- refquery@meta.data[colnames(aml_dx.obj_merged), "predicted.celltype"] 

saveRDS(refquery, file = "/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_merged_with_Lasry_et_al.rds")









