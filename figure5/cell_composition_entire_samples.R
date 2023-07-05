#### Figure 6 
library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(stringr)

"%ni%" = Negate("%in%")
options(ggrepel.max.overlaps = Inf)
output_dir <- "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/manuscript/figure6/"

all.combined.new_filtered <- readRDS("/gpfs/home/sl7424/sl7424/AML_IDH_NEW/seurat_objects/final_objects/all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_adt_processed_092622.rds")
all.combined.new_filtered$manual_annotation_v3_l3_final <- as.character(all.combined.new_filtered$manual_annotation_v3_l1_final)
all.combined.new_filtered$manual_annotation_v3_l3_final[all.combined.new_filtered$manual_annotation_v3_l1_final %in% c("HSC", "MPP", "Immature Mono")] <- "HSC_MPP"
all.combined.new_filtered$Source[all.combined.new_filtered$donor.id %in% c(paste0("BM_healthy_NYU_", 1:8), paste0("healthy_MSK_0", 1:2))] <- "BM"
all.combined.new_filtered$Source[all.combined.new_filtered$donor.id %in% c(paste0("healthy_MSK_0", 3:4))] <- "PB"

## make a data frame 
cell.numbers <- data.frame(cell_type= all.combined.new_filtered@meta.data$manual_annotation_v3_l3_final, 
                           stage = all.combined.new_filtered@meta.data$stage_simple, 
                           donor.id = all.combined.new_filtered@meta.data$donor.id,
                           Source = all.combined.new_filtered$Source,
                           therapy = all.combined.new_filtered$therapy)

cell.numbers_subset <- subset(cell.numbers, stage %in% c("healthy" ,"diagnosis","remission", "relapse")) %>% 
  subset(Source %in% "BM") %>% 
  subset(therapy %ni% c("chemo"))

cell.numbers_subset <- subset(cell.numbers_subset, cell_type %in% c("HSC_MPP", "MEP", "GMP", "CD14+ Mono"))
cell.numbers_subset$cell_type <- factor(cell.numbers_subset$cell_type, levels = c("HSC_MPP", "MEP", "GMP", "CD14+ Mono"))

## subset the stem & myeloid cell types 
df.cell.numbers <- data.frame(cell.numbers_subset) %>% 
  group_by(donor.id, cell_type, stage, therapy) %>% 
  summarise(count=n()) %>% 
  group_by(donor.id, stage) %>%
  mutate(perc=count/sum(count))


ggplot(df.cell.numbers, aes(x= stage, y= perc, fill = stage)) +
  geom_boxplot(outlier.shape=NA, alpha = 0.8) +
  geom_point(aes(x= stage, y= perc, fill = stage), size = 0.8,
             data = df.cell.numbers %>% subset(stage != "healthy")) +
  geom_point(aes(x= stage, y= perc, fill = stage), size = 0.8, 
             data = df.cell.numbers %>% subset(stage == "healthy")) +
  geom_line(aes(group = donor.id), color='gray', alpha=0.6, 
            data = df.cell.numbers %>% subset(stage != "healthy")) +
  scale_fill_manual(values = c("grey", brewer.pal(4, "Dark2")[1:4]))+
  facet_wrap(~cell_type, ncol=5)+ 
  theme_classic() + RotatedAxis()+
  xlab("") +
  ylab("Ratio")

ggsave(paste0(output_dir, "/cell_composition_boxplots_relapse_patients.pdf"), width =6, height =3, dpi = 300, device = "pdf")

