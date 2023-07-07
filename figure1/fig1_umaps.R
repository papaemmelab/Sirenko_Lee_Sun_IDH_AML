### Plot UMAPs from Figure 1: 

library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyr)
# library(EnhancedVolcano)
library(purrr)
# library(clustree)
library(SeuratDisk)
# library(msigdbr)
library(SeuratData)
source("../plotting_helpers.R")

int_obj_all = readRDS(paste0(outdir_all, "int_obj_all_CRv5_20220824_scores_GoT_numbat_20220812.rds"))


### merge umap coord 
meta_umap = as.data.frame(Embeddings(object = int_obj_all[["rpca_umap"]]))
colnames(meta_umap) = c("rpca_UMAP_1", "rpca_UMAP_2")

int_obj_all = AddMetaData(int_obj_all, metadata=meta_umap)
meta_all = int_obj_all[[]]


### Plot umap all 

meta_all0 = meta_all %>% dplyr::filter(!(donor.id == "CH-06"))

# remove weird cells 
meta_all0 = meta_all0 %>% mutate(clean = case_when((manual_annotation_v3_l1_final %in% c("B")) & (rpca_UMAP_2 > -10) ~ "remove",
                                                    !(manual_annotation_v3_l1_final %in% c("CD4+ T","CD8+ T","NK" )) & (rpca_UMAP_1 > 5) & (rpca_UMAP_2 > -5) ~ "remove",
                                                    TRUE ~ "keep"

)) %>% dplyr::filter(clean == "keep")

umap_pt_size = 0.001



### COLOR UMAP BY CELL TYPE ANNOTATION 
p_all_umap = ggplot(meta_all0,  aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + 
  geom_point( aes(color = as.factor(manual_annotation_v3_l1_final)) , size=umap_pt_size) + 
  scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
  theme_umap() +  
  ggtitle("All UMAP") + ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = 2)))
save_plot(p_all_umap, "umap_all_celltype", 8, 8)


### COLOR UMAP BY PATIENT 

shuf <- function(df){
  return(df[sample(1:dim(df)[1], dim(df)[1]),])
}

shuf_meta_all = shuf(meta_all0) # shuffle cells for plotting order 
p_all_umap = ggplot(shuf_meta_all,  aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + 
  geom_point(data = subset(shuf_meta_all, (institute_disease %in% c("NYU_healthy", "MSK_healthy"))), colour = makeTransparent("grey") , size=umap_pt_size) + 
  geom_point(data = subset(shuf_meta_all, (institute_disease %in% c("NYU_AML", "MSK_AML", "MSK_CH"))), aes(color = as.factor(donor.id)) , 
  size=umap_pt_size)+ scale_color_manual(values = patient_pal_v2, na.value = makeTransparent("grey"))   + 
  theme_umap() +  
  ggtitle("All UMAP") + ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = 2)))
save_plot(p_all_umap, "umap_all_color_patient", 8, 8)


