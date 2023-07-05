#### Figure 6 
library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(cowplot)

"%ni%" = Negate("%in%")
options(ggrepel.max.overlaps = Inf)
output_dir <- "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/manuscript/figure6/"

pal_use_man = c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
                "#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660", "#D7B5A6", "#FF9D9A" ,#"#79706E" ,"#BAB0AC",# '#0f5959', '#17a697', '#638ca6')
                "#A0CBE8", "#FABFD2" , "#B07AA1") #, '#380456', '#2a1c72')
#  '#2463a7', '#5f8fff', '#a76350', '#ce6f4a', ) #, '#d93240' "#499894",)
names(pal_use_man) = c("HSC",  "GMP" , "MPP" ,  "MEP" , "Ery", "Progenitor_DC" , "cDC", "CD14+ Mono" ,  "CD16+ Mono" ,
                       "Macrophage",  "Progenitor_B" ,  "Plasma"  , "Stroma" , 
                       "B" ,  "CD4+ T" ,"CD8+ T", "Plasmablast" , "NK", "Platelet", "pDC")

all.combined.new_filtered <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_adt_processed_092622.rds")
all.combined.new_filtered$manual_annotation_v3_l3_final <- as.character(all.combined.new_filtered$manual_annotation_v3_l1_final)
all.combined.new_filtered$manual_annotation_v3_l3_final[all.combined.new_filtered$manual_annotation_v3_l1_final %in% c("Immature Mono")] <- "MPP"

all.combined.new_filtered$Source[all.combined.new_filtered$donor.id %in% c(paste0("BM_healthy_NYU_", 1:8), paste0("healthy_MSK_0", 1:2))] <- "BM"
all.combined.new_filtered$Source[all.combined.new_filtered$donor.id %in% c(paste0("healthy_MSK_0", 3:4))] <- "PB"

all.combined.new_filtered$stage <- factor(all.combined.new_filtered$stage, 
                                          levels= c("healthy", "Dx","PR", "CR", "CR1", "CR2", "R", "R1", "R2", "P1"))

all.combined.new_filtered$stage_simple_v2 <- as.character(all.combined.new_filtered$stage_simple)
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage == "Dx" & all.combined.new_filtered$donor.id == "AML2123"] <- "diagnosis"
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage %in% c("CR1", "CR2")& all.combined.new_filtered$donor.id == "AML2123"] <- "remission"
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage %in% c("R1", "R2")& all.combined.new_filtered$donor.id == "AML2123"] <- "relapse1"
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage == "P1"& all.combined.new_filtered$donor.id == "AML2123"] <- "relapse2"

all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage == "Dx" & all.combined.new_filtered$donor.id == "AML0048"] <- "diagnosis"
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage == c("R1")& all.combined.new_filtered$donor.id == "AML0048"] <- "relapse1"
all.combined.new_filtered$stage_simple_v2[all.combined.new_filtered$stage == "R2"& all.combined.new_filtered$donor.id == "AML0048"] <- "relapse2"

all.combined.new_filtered$stage_simple_v2 <- factor(all.combined.new_filtered$stage_simple_v2, levels = c("healthy", "diagnosis", "remission",
                                                                      "relapse1", "relapse2"))
#############################################################################################################################################
### plot for AML2123 
#############################################################################################################################################

donor.id = "AML2123"
celltypes = c("HSC", "GMP", "MPP", "MEP", "Ery", "Progenitor_DC", "cDC", "CD14+ Mono", "CD16+ Mono", "Macrophage")
seurat_object <- all.combined.new_filtered[, all.combined.new_filtered$donor.id %in% c(donor.id, paste0("BM_healthy_NYU_", c(1:8)))]
seurat_object <- seurat_object[, seurat_object$manual_annotation_v3_l3_final %in% celltypes]

## make a data frame 
umap_df <- data.frame(seurat_object[[reduction]]@cell.embeddings,
                      stage = seurat_object[["stage_simple_v2"]][[1]])
umap_df$celltype <- as.character(seurat_object$manual_annotation_v3_l3_final)
umap_df$celltype[umap_df$stage == "healthy"] <- "healthy"

## save the umaps 
p <- list()
stage_plot <- c("diagnosis", "remission", "relapse1", "relapse2")
for (i in 1:length(stage_plot)){
  p[[i]] <- ggplot(umap_df, aes(RPCAUMAP_1, RPCAUMAP_2))+
    geom_point(data= umap_df %>% subset(stage == "healthy"), 
               color = "lightgrey", alpha =0.5, size = 0.2) + 
    geom_point(data= umap_df %>% subset(stage == stage_plot[i]), 
               aes(color = celltype), size = 0.2) +
    theme_classic()+
    theme(legend.position = "none") +
    scale_color_manual(values = pal_use_man, name = "celltype")+
    xlab("UMAP 1") + ylab("UMAP 2")
}

### stacked barplots for quantification 
df.cell.numbers <- data.frame(seurat_object@reductions$rpca_umap@cell.embeddings)
df.cell.numbers$celltype <- as.character(seurat_object$manual_annotation_v3_l3_final)
df.cell.numbers$stage <- seurat_object$stage_simple_v2
df.cell.numbers <- data.frame(df.cell.numbers) %>% 
  group_by(celltype, stage) %>% 
  summarise(count=n()) %>% 
  #group_by(sample) %>%
  group_by(stage) %>%
  mutate(perc=count/sum(count)) 

p[[5]] <- ggplot(df.cell.numbers, aes(x= stage, y= perc, fill = celltype)) +
  geom_col()+
  scale_fill_manual(values = pal_use_man, name = NULL)+
  theme_classic() + RotatedAxis() + xlab("") + 
  guides(fill=guide_legend(ncol=2))

## save the plots 
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], rel_widths = c(1,1,1,1,1.5),  aligh = "hv", nrow=1, ncol=5)
ggsave(paste0(output_dir, donor.id,  "_UMAP", "_" , "celltype_stage", ".pdf" ), width =18, height =3.5, dpi = 300, device = "pdf")

#############################################################################################################################################
### plot for AML0048
#############################################################################################################################################
donor.id = "AML0048"
celltypes = c("HSC", "GMP", "MPP", "MEP", "Ery", "Progenitor_DC", "cDC", "CD14+ Mono", "CD16+ Mono", "Macrophage")
seurat_object <- all.combined.new_filtered[, all.combined.new_filtered$donor.id %in% c(donor.id, paste0("BM_healthy_NYU_", c(1:8)))]
seurat_object <- seurat_object[, seurat_object$manual_annotation_v3_l3_final %in% celltypes]

## make a data frame 
umap_df <- data.frame(seurat_object[[reduction]]@cell.embeddings,
                      stage = seurat_object[["stage_simple_v2"]][[1]])
umap_df$celltype <- as.character(seurat_object$manual_annotation_v3_l3_final)
umap_df$celltype[umap_df$stage == "healthy"] <- "healthy"

## save the umaps 
p <- list()
stage_plot <- c("diagnosis", "remission", "relapse1", "relapse2")
for (i in 1:length(stage_plot)){
  p[[i]] <- ggplot(umap_df, aes(RPCAUMAP_1, RPCAUMAP_2))+
    geom_point(data= umap_df %>% subset(stage == "healthy"), 
               color = "lightgrey", alpha =0.5, size = 0.2) + 
    geom_point(data= umap_df %>% subset(stage == stage_plot[i]), 
               aes(color = celltype), size = 0.2) +
    theme_classic()+
    theme(legend.position = "none") +
    scale_color_manual(values = pal_use_man, name = "celltype")+
    xlab("UMAP 1") + ylab("UMAP 2")
}

### stacked barplots for quantification 
df.cell.numbers <- data.frame(seurat_object@reductions$rpca_umap@cell.embeddings)
df.cell.numbers$celltype <- as.character(seurat_object$manual_annotation_v3_l3_final)
df.cell.numbers$stage <- seurat_object$stage_simple_v2
df.cell.numbers <- data.frame(df.cell.numbers) %>% 
  group_by(celltype, stage) %>% 
  summarise(count=n()) %>% 
  #group_by(sample) %>%
  group_by(stage) %>%
  mutate(perc=count/sum(count)) 

p[[5]] <- ggplot(df.cell.numbers, aes(x= stage, y= perc, fill = celltype)) +
  geom_col()+
  scale_fill_manual(values = pal_use_man, name = NULL)+
  theme_classic() + RotatedAxis() + xlab("") + 
  guides(fill=guide_legend(ncol=2))

## save the plots 
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], rel_widths = c(1,1,1,1,1.5),  aligh = "hv", nrow=1, ncol=5)
ggsave(paste0(output_dir, donor.id,  "_UMAP", "_" , "celltype_stage", ".pdf" ), width =18, height =3.5, dpi = 300, device = "pdf")



