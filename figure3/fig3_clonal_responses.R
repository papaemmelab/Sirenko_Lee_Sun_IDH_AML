library(Seurat)
library(dplyr)
library(tidyverse)
library(tidyr)
source("../plotting_helpers.R")



### Pre - Post - treatment figure 
df_plot_pre = int_obj_all[[]] %>% dplyr::filter(Timepoint == "T1", institute_disease == "MSK_AML") %>% mutate(group = "AML")
df_plot_post = int_obj_all[[]] %>% dplyr::filter(Timepoint == "T2") %>% mutate(group = "AML")
df_plot_norm = int_obj_all[[]]  %>% dplyr::filter( institute_disease %in% c("MSK_healthy", "NYU_healthy")) %>% mutate(group = "Normal")
df_plot_pre = rbind(df_plot_pre, df_plot_norm)
df_plot_post = rbind(df_plot_post, df_plot_norm)
patients_plot_MSK = c( "IDH1i-01", "IDH1i-02", "IDH1i-03" ,"IDH2i-01", "IDH2i-02", "IDH2i-03")
for (i in 1:length(names_dict)){patients_plot_MSK <- replace(patients_plot_MSK, patients_plot_MSK == names(names_dict[i]), as.character(names_dict[i]))}
for (i in 1:length(names_dict)){df_plot_pre$donor.id <- replace(df_plot_pre$donor.id, df_plot_pre$donor.id == names(names_dict[i]), as.character(names_dict[i]))}
for (i in 1:length(names_dict)){df_plot_post$donor.id <- replace(df_plot_post$donor.id, df_plot_post$donor.id == names(names_dict[i]), as.character(names_dict[i]))}


df_plot_pre = df_plot_pre %>% mutate(mutant_any = case_when( (IDH2_R140 == "MT" | IDH2_R172 == "MT"  | NRAS_G12 == "MT" | NPM1_W288 == "MT" | SRSF2_P95 == "MT" |  IDH1_R132 == "MT" | clone %in% c("gain_10","gain_14", "gain_1q_8", "gain_6_10","gain_8" )) ~ "MT",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT" | clone %in% "normal" ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_any = factor(mutant_any, levels = c("MT", "WT", "NA")))
df_plot_post = df_plot_post %>% mutate(mutant_any = case_when( (IDH2_R140 == "MT" | IDH2_R172 == "MT"  | NRAS_G12 == "MT" | NPM1_W288 == "MT" | SRSF2_P95 == "MT" |  IDH1_R132 == "MT" | clone %in% c("gain_10","gain_14", "gain_1q_8", "gain_6_10","gain_8" )) ~ "MT",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT" | clone %in% "normal" ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_any = factor(mutant_any, levels = c("MT", "WT", "NA")))


df_plot_pre = df_plot_pre %>% mutate(mutant_any_GoT = case_when( (IDH2_R140 == "MT" | IDH2_R172 == "MT"  | NRAS_G12 == "MT" | NPM1_W288 == "MT" | SRSF2_P95 == "MT" |  IDH1_R132 == "MT" ) ~ "MT",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT"  ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_any = factor(mutant_any, levels = c("MT", "WT", "NA")))
df_plot_post = df_plot_post %>% mutate(mutant_any_GoT = case_when( (IDH2_R140 == "MT" | IDH2_R172 == "MT"  | NRAS_G12 == "MT" | NPM1_W288 == "MT" | SRSF2_P95 == "MT" |  IDH1_R132 == "MT" ) ~ "MT",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT"  ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_any = factor(mutant_any, levels = c("MT", "WT", "NA")))

table(df_plot_pre$mutant_any_GoT) # (14824+3070) / 77364 = 0.2312962
table(df_plot_post$mutant_any_GoT)

## color by mutation 
df_plot_pre = df_plot_pre %>% mutate(mutant_specific = case_when(NRAS_G12 == "MT" ~ "NRAS", 
                                                              NPM1_W288 == "MT" ~ "NPM1",
                                                              SRSF2_P95 == "MT" ~ "SRSF2", 
                                                              (IDH2_R140 == "MT" | IDH2_R172 == "MT"  ) ~ "IDH2", 
                                                              IDH1_R132 == "MT" ~ "IDH1", 
                                                              clone %in% c("gain_10", "gain_6_10") ~ "gain_6_10", 
                                                              clone %in% c("gain_14") ~ "gain_14", 
                                                              clone %in% c("gain_1q_8") ~ "gain_1q_8", 
                                                              clone %in% c("gain_8" ) ~ "gain_8",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT" ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_specific = factor(mutant_specific, levels = c("NRAS", "SRSF2", "NPM1", "IDH2", "IDH1", "gain_6_10", "gain_14", "gain_1q_8", "gain_8","WT", "NA")))

df_plot_post = df_plot_post %>% mutate(mutant_specific = case_when(NRAS_G12 == "MT" ~ "NRAS", 
                                                              NPM1_W288 == "MT" ~ "NPM1",
                                                              SRSF2_P95 == "MT" ~ "SRSF2", 
                                                              (IDH2_R140 == "MT" | IDH2_R172 == "MT"  ) ~ "IDH2", 
                                                              IDH1_R132 == "MT" ~ "IDH1", 
                                                              clone %in% c("gain_10", "gain_6_10") ~ "gain_6_10", 
                                                              clone %in% c("gain_14") ~ "gain_14", 
                                                              clone %in% c("gain_1q_8") ~ "gain_1q_8", 
                                                              clone %in% c("gain_8" ) ~ "gain_8",
                                                     (IDH2_R140 == "WT" | IDH2_R172 == "WT"  | NRAS_G12 == "WT" | NPM1_W288 == "WT" | SRSF2_P95 == "WT"| IDH1_R132 == "WT" ) ~ "WT",
                                                    (TRUE ~ "NA")) ) %>% mutate(mutant_specific = factor(mutant_specific, levels = c("NRAS", "SRSF2", "NPM1", "IDH2", "IDH1", "gain_6_10", "gain_14", "gain_1q_8", "gain_8","WT", "NA")))

# pal_mut_specific = c("#11688b", "#0c4f4e", "#9d442a", "#5c315f", "#b64c63", "#008080", "orange", "pink", "yellow", "beige", "grey")
# names(pal_mut_specific) = c("NRAS", "SRSF2", "NPM1", "IDH2", "IDH1", "gain_6_10", "gain_14", "gain_1q_8", "gain_8","WT", "NA")

pal_mut_specific = c("#11688b", "#0c4f4e", "#AA7A38", "#5c315f", "#b64c63", "#008080", "#FFBE7D", "#D4A6C8", "#B07AA1", "#F1CE63", "grey")
names(pal_mut_specific) = c("NRAS", "SRSF2", "NPM1", "IDH2", "IDH1", "gain_6_10", "gain_14", "gain_1q_8", "gain_8","WT", "NA")


co_mut_colors = c("#b64c63", "#5c315f", "#AA7A38", "#11688b", "#0c4f4e")
names(co_mut_colors) = c("IDH1" , "IDH2",  "NPM1" , "NRAS" , "SRSF2")

## filter for LD celltypes 
df_plot_pre$cell_type = df_plot_pre$manual_annotation_v3_l1_final
df_plot_post$cell_type = df_plot_post$manual_annotation_v3_l1_final

df_plot_pre_LD = df_plot_pre %>% dplyr::filter(cell_type %in% c("CD14+ Mono","CD16+ Mono", "MEP","Immature Mono", "cDC",  "GMP" , "MPP" ,"Ery", "HSC", "Progenitor_DC"))
df_plot_post_LD = df_plot_post %>% dplyr::filter(cell_type %in% c("CD14+ Mono","CD16+ Mono", "MEP","Immature Mono", "cDC",  "GMP" , "MPP" ,"Ery", "HSC", "Progenitor_DC"))


# remove weird cells 
df_plot_pre_LD = df_plot_pre_LD %>% mutate(clean = case_when((cell_type %in% c("B")) & (rpca_UMAP_2 > -10) ~ "remove",
                                                    !(cell_type %in% c("CD4+ T","CD8+ T","NK" )) & (rpca_UMAP_1 > 5) & (rpca_UMAP_2 > -5) ~ "remove",
                                                    TRUE ~ "keep"

)) %>% dplyr::filter(clean == "keep")
df_plot_post_LD = df_plot_post_LD %>% mutate(clean = case_when((cell_type %in% c("B")) & (rpca_UMAP_2 > -10) ~ "remove",
                                                    !(cell_type %in% c("CD4+ T","CD8+ T","NK" )) & (rpca_UMAP_1 > 5) & (rpca_UMAP_2 > -5) ~ "remove",
                                                    TRUE ~ "keep"

)) %>% dplyr::filter(clean == "keep")


p_umap_post = list()
point_size = 0.001
figh_in = 20
figh_in_MSK = 1.5*6

for (i in 1:length(patients_plot_MSK)){
  p_umap_post[[i]] = ggplot(df_plot_post_LD,  mapping = aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
  geom_point(data = subset(df_plot_post_LD, mutant_specific %in% "NA"),  fill = makeTransparent("grey") ,size=point_size,colour=makeTransparent("grey")) + 
  geom_point(subset(subset(df_plot_post_LD, donor.id %in% patients_plot_MSK[i]), mutant_specific %in% c( "IDH2", "IDH1","gain_14", "gain_1q_8", "gain_8","WT")), mapping =aes( colour = mutant_specific) ,size=point_size)+ #, colour = manual_celltype.l3
  geom_point(subset(subset(df_plot_post_LD, donor.id %in% patients_plot_MSK[i]), mutant_specific %in% c("NRAS", "SRSF2", "NPM1", "gain_6_10", "gain_14")), mapping =aes( colour = mutant_specific) ,size=point_size)+ #, colour = manual_celltype.l3
  scale_fill_manual(values = pal_use_man, na.value = makeTransparent("grey")) + scale_color_manual(values = pal_mut_specific, na.value = makeTransparent("grey")) + 
   xlab("") + ylab ("") + ggtitle(patients_plot_MSK[i]) +  theme(legend.position='none') + labs(colour="Mutant") +guides(color = guide_legend(override.aes = list(size = leg_size)))+ 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  }

cow_celltype_post = cowplot::plot_grid(plotlist = p_umap_post, ncol = 1)
ggsave(plot = cow_celltype_post, width = 1.5, height = figh_in_MSK, filename = paste0(outdir_all, "figures_integrated_object/umap_T2_mutant_specific_nolegend_fig3_LD.png"))
ggsave(plot = cow_celltype_post, width = 1.5, height = figh_in_MSK, filename = paste0(outdir_all, "figures_integrated_object/umap_T2_mutant_specific_nolegend_fig3_LD.pdf"))

p_umap_pre = list()
for (i in 1:length(patients_plot_MSK)){
  p_umap_pre[[i]] = ggplot(df_plot_pre_LD,  mapping = aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
  geom_point(data = subset(df_plot_pre_LD, mutant_specific %in% "NA"),  fill = makeTransparent("grey") ,size=point_size,colour=makeTransparent("grey")) + 
  geom_point(subset(subset(df_plot_pre_LD, donor.id %in% patients_plot_MSK[i]), mutant_specific %in% c( "IDH2", "IDH1","gain_14", "gain_1q_8", "gain_8","WT")), mapping =aes( colour = mutant_specific) ,size=point_size)+ #, colour = manual_celltype.l3
  geom_point(subset(subset(df_plot_pre_LD, donor.id %in% patients_plot_MSK[i]), mutant_specific %in% c("NRAS", "SRSF2", "NPM1", "gain_6_10", "gain_14")), mapping =aes( colour = mutant_specific) ,size=point_size)+ #, colour = manual_celltype.l3
  scale_fill_manual(values = pal_use_man, na.value = makeTransparent("grey")) + scale_color_manual(values = pal_mut_specific, na.value = makeTransparent("grey")) + 
  xlab("") + ylab ("") + ggtitle(patients_plot_MSK[i]) + theme(legend.position='none') + labs(colour="Mutant") +# guides(color = guide_legend(override.aes = list(size = leg_size)))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  }

cow_celltype_pre = cowplot::plot_grid(plotlist = p_umap_pre, ncol = 1)
ggsave(plot = cow_celltype_pre, width = 1.5, height = figh_in_MSK, filename = paste0(outdir_all, "figures_integrated_object/umap_T1_mutant_specific_no_legend_fig3_LD.png"))
ggsave(plot = cow_celltype_pre, width = 1.5, height = figh_in_MSK, filename = paste0(outdir_all, "figures_integrated_object/umap_T1_mutant_specific_no_legend_fig3_LD.pdf"))



## stacked bar plot pre/post treatment mutant_specific 
meta_data = rbind(df_plot_pre_LD, df_plot_post_LD)
meta_data_stack = meta_data %>% dplyr::filter(Patient %in% c("IDH1i-01" ,"IDH1i-02" ,"IDH1i-03", "IDH2i-01","IDH2i-02","IDH2i-03"), Timepoint %in% c("T1", "T2"), cell_type %in%c("CD14+ Mono","CD16+ Mono", "MEP","Immature Mono", "cDC",  "GMP" , "MPP" ,"Ery", "HSC" )) # , "Progenitor_DC"  ))


meta_data_stack = meta_data_stack %>% mutate(cell_type_broad_simple = case_when(cell_type %in% c("HSC", "MPP", "Immature Mono") ~ "HSC/MPP",
                                                                          cell_type %in% c("GMP") ~ "GMP", 
                                                                          cell_type %in% c("CD14+ Mono" , "CD16+ Mono", "cDC") ~ "Mono/DC",
                                                                          cell_type %in% c("MEP", "Ery") ~ "MEP/Ery", 
                                                                          # cell_type %in% c("Immature Mono") ~ "Imm Mono"
                                                                          )) %>% mutate(cell_type_broad_simple = factor(cell_type_broad_simple, levels = c( "MEP/Ery" ,  "Mono/DC", "Imm Mono","GMP","HSC/MPP" ))) %>% 
                                                                          mutate(mutant_specific <- factor(mutant_specific, levels=c("NA", "WT", "gain_14" ,  "gain_1q_8" ,"gain_8",   "IDH2","IDH1", "NRAS","SRSF2","NPM1",  "gain_6_10")))

## get values for text 
pt = c("AML-06" ) #,"AML-02" ,"AML-03", "AML-04","AML-05","AML-06")
  df2 <- meta_data_stack %>% dplyr::filter(donor.id == pt & Timepoint == "T1") %>% 
  group_by(cell_type_broad_simple, Timepoint, mutant_any) %>%
  summarise(n=n()) 
  

## plot stacked bar plots for Figure 3
p_stack_list = list()
p_stack_list_noleg = list()
p_stack_list_leg = list()
for (pt in c("AML-01" ,"AML-02" ,"AML-03", "AML-04","AML-05","AML-06")){

# # calc n observations 
  df2 <- meta_data_stack %>% dplyr::filter(donor.id == pt) %>% 
  group_by(cell_type_broad_simple, Timepoint) %>%
  summarise(n=n()) %>% dplyr::mutate(xpos = case_when(Timepoint == "T1" ~ -0.5, Timepoint == "T2" ~ 0.5)) 
  
p_stack_list[[pt]] = meta_data_stack %>% dplyr::filter(donor.id == pt) %>% 
  group_by(cell_type_broad_simple, mutant_specific) %>% 
  summarize(T1 = -sum(Timepoint == "T1"), T2 = sum(Timepoint == "T2")) %>%
  ggplot() + 
    geom_col(aes(y = cell_type_broad_simple, x = T1, fill = factor(mutant_specific, levels=c("NA", "WT", "gain_14" ,  "gain_1q_8" ,"gain_8",   "IDH2","IDH1", "NRAS","SRSF2","NPM1",  "gain_6_10"))), position = "fill") +
    geom_col(aes(y = cell_type_broad_simple, x = T2, fill = factor(mutant_specific, levels=c("NA", "WT", "gain_14" ,  "gain_1q_8" ,"gain_8",   "IDH2","IDH1", "NRAS","SRSF2","NPM1",  "gain_6_10"))), position = "fill") + 
    geom_hline(aes(yintercept = 0)) + 
    labs(x = paste("T1", "T2", sep = paste(rep(" ", 6), collapse = " ")), y = "", title = pt) + 
    geom_text(data = df2, aes(y = cell_type_broad_simple, x = xpos, label = paste0("n=",n))) +
    # geom_text(aes(label=after_stat('count')), #stat='count', 
    # nudge_y=0.125) + 
    # geom_text(stat = "count", position = position_fill(.5)) + 
    
    theme_classic() + scale_fill_manual(values = pal_mut_specific, limits = force) + 
    geom_vline(xintercept = 0, color = "black", size=0.5)  +
    scale_x_continuous(breaks=c(-1, -0.5, 0, 0.5, 1), labels=c("1", "0.5", " ", "0.5", "1"))  
    # geom_text(aes(label=after_stat(count))) 
    # geom_label(stat = 'count', aes(y = manual_celltype.l3, label = after_stat(count)))
    
    # geom_label(stat = 'count', aes(label = ..count..), 
    #          vjust = -0.1,
    #          show.legend = FALSE) 
    # geom_label(aes(label = after_stat(count)), stat = "count")

p_stack_list_noleg[[pt]] = p_stack_list[[pt]] + theme(legend.position="none") + theme(legend.title=element_blank())
p_stack_list_leg[[pt]] = p_stack_list[[pt]] + theme(legend.position="left") + theme(legend.title=element_blank())

}
p_grid_mut <- cowplot::plot_grid(plotlist = p_stack_list_noleg , ncol = 1, nrow = 6)
ggsave(plot = cowplot::plot_grid(p_grid_mut), width = 3, height = 10, filename = paste0(outdir_all, "figures_integrated_object/stacked_barplot_inverted_mutant_specific_cell_type_broad_simple_fig3.png"))
ggsave(plot = cowplot::plot_grid(p_grid_mut), width = 3, height = 10, filename = paste0(outdir_all, "figures_integrated_object/stacked_barplot_inverted_mutant_specific_cell_type_broad_simple_fig3.pdf"))

p_grid_mut <- cowplot::plot_grid(plotlist = p_stack_list_leg , ncol = 1, nrow = 6)
ggsave(plot = cowplot::plot_grid(p_grid_mut), width = 3, height = 10, filename = paste0(outdir_all, "figures_integrated_object/stacked_barplot_inverted_mutant_specific_cell_type_broad_simple_fig3_legend.png"))
ggsave(plot = cowplot::plot_grid(p_grid_mut), width = 3, height = 10, filename = paste0(outdir_all, "figures_integrated_object/stacked_barplot_inverted_mutant_specific_cell_type_broad_simple_fig3_legend.pdf"))




### boxplot of fraction immature / mature myeloid cells pre/post treatment MSK cohort 

df_mye_ratio_total = meta_data_stack %>% mutate(myeloid = case_when(manual_annotation_v3_l1_final %in% c("CD14+ Mono","CD16+ Mono", "cDC" ,"Macrophage"  ) ~ "mature",
                                                            manual_annotation_v3_l1_final %in% c("HSC","MPP","GMP","MEP","Immature Mono", "Progenitor_DC", "Ery" ) ~ "immature")) %>% 
                                                             dplyr::count(donor.id, myeloid, Timepoint) %>% group_by (donor.id, Timepoint) %>% summarise(mye_ratio_total = n[myeloid=="immature"]/sum(n))
df_mye_ratio_total = df_mye_ratio_total %>% mutate(donor.id = factor(donor.id, levels = rev(c("AML-06", "AML-05", "AML-04", "AML-03", "AML-02", "AML-01"))))

p_box_mye = df_mye_ratio_total %>% 
                                            ggplot(aes(x = Timepoint, y = mye_ratio_total)) + geom_boxplot(fill = "#aa6f73", outlier.shape = NA, width = 0.5, alpha = 0.3) +
                                            #scale_fill_manual(values = c("#aa6f73", "#f6e0b5")) + 
                                            geom_line(aes(group = donor.id, color = donor.id), size = 0.3) +  geom_point(aes(color = donor.id), size=0.8, alpha=0.9) +
                                            #stat_compare_means(paired =TRUE) + 
                                            scale_colour_brewer(palette = "Dark2", name = "Patient") + ylim(0,1) + 
                                            theme_classic() + scale_x_discrete(breaks=c( "T1", "T2"),labels=c("Diagnosis", "Post-Treatment")) + ylab("Fraction Immature / Total Myeloid") + xlab("") + 
                                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
save_plot_simple(p_box_mye, "boxplot_myeloid_ratio_MSK_cohort_dx_post_stat", 2.5, 4.5)                                          


t1_vals = df_mye_ratio_total %>% dplyr::filter(Timepoint == "T1") 
t2_vals = df_mye_ratio_total %>% dplyr::filter(Timepoint == "T2") 

t.test(t1_vals$mye_ratio_total, t2_vals$mye_ratio_total, paired = TRUE, alternative = "greater")


## get stats for how many cells have CNAs 

df = int_obj_all[[]] %>% dplyr::filter(donor.id == "IDH1i-02" ,  stage_simple == "diagnosis") 
table(df$clone)

df = int_obj_all[[]] %>% dplyr::filter(donor.id == "IDH2i-02" ,  stage_simple == "diagnosis") 
table(df$clone)


df = int_obj_all[[]] %>% dplyr::filter(donor.id == "IDH2i-01" ,  stage_simple == "diagnosis") 
table(df$clone)

df = int_obj_all[[]] %>% dplyr::filter(donor.id == "IDH2i-01" ,  stage_simple == "diagnosis") 





## GENE EXPRESSION VS GoT EFFICIENCY 

# calc expression
got_cohort = subset(int_obj_all, institute_disease %in% c( "MSK_CH"  ,   "MSK_AML"))
exp_df = AverageExpression(got_cohort,  features = c("IDH1", "IDH2", "DNMT3A", "SRSF2", "NRAS", "NPM1"), assay = "RNA", group.by = "institute")
exp_df0 = as.data.frame(exp_df$RNA)
exp_df0 = exp_df0 %>% tibble::rownames_to_column("gene")
colnames(exp_df0) = c("gene", "rna_exp")
# calc got efficiency 
got_res = meta_all0 %>% dplyr::select(donor.id, IDH1_R132, IDH2_R140, IDH2_R172, NPM1_W288, NRAS_G12, SRSF2_P95)
SRSF2_eff = got_res %>% dplyr::filter(donor.id %in% c("AML-04", "AML-05", "AML-06", "CH-02")) %>% select(SRSF2_P95) %>% table(useNA = "always")
NRAS_eff = got_res %>% dplyr::filter(donor.id %in% c("AML-02")) %>% select(NRAS_G12) %>% table(useNA = "always")
NPM1_eff = got_res %>% dplyr::filter(donor.id %in% c("AML-01", "AML-03")) %>% select(NPM1_W288) %>% table(useNA = "always")
IDH2_R140_eff = got_res %>% dplyr::filter(donor.id %in% c( "AML-05", "AML-06", "CH-01", "CH-02")) %>% select(IDH2_R140) %>% table(useNA = "always")
IDH2_R172_eff = got_res %>% dplyr::filter(donor.id %in% c( "AML-04")) %>% select(IDH2_R172) %>% table(useNA = "always")
IDH1_eff = got_res %>% dplyr::filter(donor.id %in% c( "AML-01", "AML-02", "AML-03")) %>% select(IDH1_R132) %>% table(useNA = "always")

got_eff = bind_rows(NRAS_eff, SRSF2_eff, NPM1_eff, IDH2_R140_eff,IDH2_R172_eff, IDH1_eff)
got_eff$gene = c("NRAS", "SRSF2", "NPM1", "IDH2_R140", "IDH2_R172", "IDH1")
colnames(got_eff) = c("MT", "WT", "none", "gene")


got_eff = got_eff %>% mutate(fraction_genotyped = (MT + WT) / (MT + WT + none))

got_eff = got_eff %>% full_join(exp_df0)
p_got = ggplot(got_eff, aes(rna_exp, fraction_genotyped, color = gene, label = gene)) + geom_point() + theme_classic() + scale_color_manual(values = co_mut_colors) + geom_text(vjust=-1) + ylab("Fraction Genotyped") + xlab("RNA Expression")
save_plot(p_got, "got_efficiency", 8, 8)

# read in distance file 
got_target_dist = read.csv(paste0(outdir_all, "got_distance_to_5prime.csv"), sep = ",", header = TRUE)
got_eff = got_eff %>% full_join(got_target_dist)
got_eff_df = as.data.frame(got_eff)
got_eff_df["rna_exp"][got_eff_df["gene"] == "IDH2_R140"] <- got_eff_df["rna_exp"][got_eff_df["gene"] == "IDH2"] 
got_eff_df["rna_exp"][got_eff_df["gene"] == "IDH2_R172"] <- got_eff_df["rna_exp"][got_eff_df["gene"] == "IDH2"] 

got_eff_df$mutation_label = gsub("_", " p.", got_eff_df$mutation )
got_eff_df = got_eff_df %>% mutate(point_label = paste0(mutation_label, "\n(", round(rna_exp, digits = 2), ")"))
p_got = got_eff_df %>% dplyr::filter(gene %in% c("NRAS", "NPM1", "SRSF2", "IDH2_R140", "IDH2_R172", "IDH1")) %>% ggplot( aes( distance, fraction_genotyped,  color = gene_name)) + geom_point(aes(size = rna_exp)) + theme_classic() + scale_color_manual(values = co_mut_colors) + 
ggrepel::geom_text_repel(aes(label = point_label), point.padding = 1) +
  scale_size(name   = "Avg RNA Exp",
             breaks = c(0.2, 1, 9),
             labels = c("0.1", "1", "10"))+ 
ylab("Fraction Cells Genotyped") + xlab("Distance from 5' End (bp)") + ylim(0, 0.4) + xlim(0, 1100)
save_plot(p_got, "got_efficiency_dist", 6, 5.5)

