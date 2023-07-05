library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(ggpubr)
library(UCell)

"%ni%" = Negate("%in%")
options(ggrepel.max.overlaps = Inf)

output_dir = "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/integrated/MSK_collaboration/diagnosis/"
mut_pal_use = c("#eaac8b", "#355070", "#6d587a", "#b56576")
names(mut_pal_use) = c("Control", "IDH", "TET2", "Others")

AML_merged <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_merged_with_Lasry_et_al.rds")

DimPlot(AML_merged, group.by = "manual_annotation_v3_l3_final" , split.by = "id",
        reduction = "rpca_umap", label = 2.5, repel = T, raster = T, pt.size = 0.5)
ggsave(paste0(output_dir, "/", "IDHm_vs_others" , "_celltype_umap.pdf"),  width =12, height =5, dpi = 300, device = "pdf")

AML_merged@meta.data[AML_merged$id == "query", "stage_simple"] <- "diagnosis"

AML_merged$celltype <- AML_merged$manual_annotation_v3_l3_final
AML_merged$celltype[AML_merged$manual_annotation_v3_l3_final%in% c("CD14+ Mono", "CD16+ Mono", "cDC")] <- "Mature_Mono_DC"
AML_merged$celltype[AML_merged$manual_annotation_v3_l3_final%in% c("HSC", "MPP")] <- "HSC_MPP"

data_df_merged <- AML_merged@meta.data %>% subset(stage_simple %in% c("healthy", "diagnosis"))
data_df_merged$cell_barcodes <- rownames(data_df_merged)
colnames(data_df_merged) <- str_remove_all(colnames(data_df_merged), pattern = "HALLMARK_|_UCell")


### plot the gene module scores 
data_df_merged_long <- reshape2::melt(data_df_merged, id.vars = c("celltype", "stage_simple", 
                                                                  "donor.id", "IDH"))

data_df_merged_long_subset <- data_df_merged_long %>% subset(variable %in% c("HYPOXIA", "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE", 
                                                                             "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
                                                                             "G2M_CHECKPOINT", "MYC_TARGETS_V2", "E2F_TARGETS", "APOPTOSIS", 
                                                                             "FATTY_ACID_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "P53_PATHWAY"
))

my_comparisons <- list( c("Control", "IDH"), c("IDH", "Other"))


df_plot <- data_df_merged_long_subset %>% subset(celltype == "Mature_Mono_DC")
ggplot(df_plot, aes(x=IDH, y= value, fill = IDH))+
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA)+
  theme_classic()+ 
  facet_wrap(~variable, ncol = 6, scales = "free") + 
  theme(strip.text = element_text(size=6),
        axis.text=element_text(size=6))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2) +
  scale_fill_manual(values = mut_pal_use, name = "Mutations")+
  ylab("Score") + xlab("") 
ggsave(paste0(output_dir, "/", "selected_genesets_score_", "Mature_Mono_DC", "_mutations_revised.pdf"),  width =11, height =4, dpi = 300, device = "pdf")



#########################################################################################################################################################
### plot the cell composition 
df.cell.numbers <- data.frame(data_df_merged) %>% 
  subset(celltype %in% c("HSC_MPP", "MEP" ,"GMP", "Mature_Mono_DC"))%>% 
  group_by(donor.id, celltype, IDH) %>% 
  summarise(count=n()) %>% 
  group_by(donor.id) %>% 
  mutate(perc=count/sum(count)) 

df.cell.numbers_subset <- subset(df.cell.numbers, donor.id != "AML1529")

df.cell.numbers_subset$celltype <- factor(df.cell.numbers_subset$celltype, 
                                          levels = c("HSC_MPP", "Mature_Mono_DC", "GMP", "MEP"))

ggplot(df.cell.numbers_subset, aes(x = IDH, y=perc, fill= IDH))+
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA)+
  theme_classic()+ 
  facet_wrap(~celltype, ncol = 4) + 
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2.5) +
  scale_fill_manual(values = mut_pal_use, name = "Mutations")+
  ylab("Ratio") + xlab("") 

ggsave(paste0(output_dir, "/", "IDHm_vs_others" , "_frequency_diagnosis_revised.pdf"),  width =9, height =3, dpi = 300, device = "pdf")





