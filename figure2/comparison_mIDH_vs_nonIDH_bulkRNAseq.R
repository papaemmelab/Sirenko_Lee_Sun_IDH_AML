#######################################################################################3
### Aliance data - deconvolution 
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

data.dir = "/gpfs/home/sl7424/sl7424/AML_IDH_NEW/Alliance_bulkRNA/Alliance_bulk_AML"
output_dir = "/gpfs/home/sl7424/sl7424/AML_IDH_NEW/Alliance_bulkRNA/results/"

# mut_pal_use = c("#355070", "#6d587a", "#b56576", "#eaac8b")
# names(mut_pal_use) = c("IDH", "TET2", "Others", "Control")
mut_pal_use = c("#355070", "#b56576", "#eaac8b")
names(mut_pal_use) = c("IDH", "Others", "Control")

alliance_meta <- read.csv(paste0(data.dir, "/mutations_for_RNAseqSamples.csv"))
alliance_deconv <- readRDS(paste0(data.dir, "/allCounts_DWLS_Alliance.rds"))
alliance_deconv[alliance_deconv<0] <- 0

df_deconv <- data.frame(HSC_MPP = rowSums(t(alliance_deconv[c("HSC", "MPP", "Immature_Mono"),])),
                        Mature_Mono_DC = rowSums(t(alliance_deconv[c("CD14pos_Mono", "CD16pos_Mono", "cDC", "Inflammatory_Mono"),])),
                        GMP = as.numeric(t(alliance_deconv[c("GMP"),])),
                        MEP = as.numeric(t(alliance_deconv[c("MEP"),])),
                        sample_id = colnames(alliance_deconv))

df_deconv_IDH <- df_deconv %>% left_join(alliance_meta[,c("sample_id", "IDH1", "IDH2", "NPM1", "SRSF2", "TET2")], by = "sample_id")
df_deconv_IDH_long <- reshape2::melt(df_deconv_IDH, id.vars = c("IDH1", "IDH2", "sample_id",  "NPM1", "SRSF2", "TET2"))

# df_deconv_IDH_long <- df_deconv_IDH_long %>% mutate(mutation = ifelse(IDH1==1 | IDH2==1, "IDH", ifelse(TET2 == 1, "TET2",  "Others")))
# df_deconv_IDH_long$mutation <- factor(df_deconv_IDH_long$mutation, levels = c( "IDH", "TET2", "Others"))
# 
# # my_comparisons <- list(c("IDH", "Others"), c("IDH", "TET2"))
# 

df_deconv_IDH_long <- df_deconv_IDH_long %>% mutate(mutation = ifelse(IDH1==1 | IDH2==1, "IDH", "Others"))
df_deconv_IDH_long$mutation <- factor(df_deconv_IDH_long$mutation, levels = c( "IDH", "Others"))
my_comparisons <- list(c("IDH", "Others"))


ggplot(df_deconv_IDH_long, aes(x = mutation, y=value, fill= mutation))+
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill = "white", outlier.shape = NA)+
  theme_classic()+ 
  facet_wrap(~variable, ncol = 6) + 
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", size = 2.5) +
  scale_fill_manual(values = mut_pal_use, name = "Mutations")+
  ylab("Ratio") + xlab("") 

ggsave(paste0(output_dir, "/", "IDHm_vs_others" , "_Alliance.pdf"),  width =8, height =3, dpi = 300, device = "pdf")


