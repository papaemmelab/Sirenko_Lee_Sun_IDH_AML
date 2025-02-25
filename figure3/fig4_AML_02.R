## Figure 4 

source("../plotting_helpers.R")

## Figure params and palettes 
NRAS_pal = c("#11688b", "#c28c09", makeTransparent("grey")) # WT #F1CE63
names(NRAS_pal)  = c("MT", "WT", "NA")

IDH1_pal = c("#b64c63", "#c28c09", makeTransparent("grey"))# WT #F1CE63
names(IDH1_pal) = c("MT", "WT", "NA")
IDH_NRAS_pal = c("#b64c63", "#11688b")

clone_pal = c("#b64c63", "#008080", "steelblue3", "beige")
names(clone_pal) = c("gain_1q_8", "gain_6_10", "gain_10", "normal")

pal_mut_specific = c("#11688b", "#0c4f4e", "#AA7A38", "#5c315f", "#b64c63", "#008080", "#FFBE7D", "#D4A6C8", "#B07AA1", "#F1CE63", "grey")
names(pal_mut_specific) = c("NRAS", "SRSF2", "NPM1", "IDH2", "IDH1", "gain_6_10", "gain_14", "gain_1q_8", "gain_8","WT", "NA")

pal_clones_manual = c("#008080", "#11688b", "#b64c63" )
names(pal_clones_manual) = c("gain_6_10_clone"  ,    "NRAS_clone",   "parent_clone" )

umap_pt_size = 0.2

## Read in seruat object 

IDH1i_02_obj = readRDS("all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_081722_subset_IDH1i-02.rds")
leg_size = 3
meta_IDH1i_02 = IDH1i_02_obj[[]]
nrow = meta_IDH1i_02 %>% mutate(cell_type_simple = case_when(manual_annotation_v3_l1_final %in% c("Immature Mono") ~ "MPP", 
                                                                      TRUE ~ as.character(manual_annotation_v3_l1_final)))


## Plot UMAPs colored by cell type and genotype

p_all_umap = ggplot(meta_IDH1i_02,  aes(x=IDH1i_02_UMAP_1, y=IDH1i_02_UMAP_2) ) + theme_minimal() + 
    geom_point( aes(color = cell_type_simple) , size=umap_pt_size)+ scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    theme_umap() + ggtitle("IDH1i-02 UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH1i-02_T1", 4, 4)

p_all_umap = ggplot(meta_IDH1i_02, aes(x=IDH1i_02_UMAP_1, y=IDH1i_02_UMAP_2, color = NRAS_G12) ) + 
    geom_point(data = subset(meta_IDH1i_02, (is.na(NRAS_G12))),  size=umap_pt_size) +  scale_color_manual(values = NRAS_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH1i_02, (NRAS_G12 == "WT")), size=umap_pt_size) +  scale_color_manual(values = NRAS_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH1i_02, (NRAS_G12 == "MT")), size=umap_pt_size) +  scale_color_manual(values = NRAS_pal, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("NRAS G12")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH1i-02-T1_CRv5_color_NRAS_G12_2", 4, 4)

p_all_umap = ggplot(meta_IDH1i_02, aes(x=IDH1i_02_UMAP_1, y=IDH1i_02_UMAP_2, color = IDH1_R132) ) + 
    geom_point(data = subset(meta_IDH1i_02, (is.na(IDH1_R132))),  size=umap_pt_size) +  scale_color_manual(values = IDH1_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH1i_02, (IDH1_R132 == "WT")), size=umap_pt_size) +  scale_color_manual(values = IDH1_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH1i_02, (IDH1_R132 == "MT")), size=umap_pt_size) +  scale_color_manual(values = IDH1_pal, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("IDH1 R132")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH1i-02-T1_CRv5_color_IDH1_R132", 4, 4)

p_all_umap = ggplot(meta_IDH1i_02, aes(x=IDH1i_02_UMAP_1, y=IDH1i_02_UMAP_2, color = clone) ) + 
    geom_point(  size=umap_pt_size) +  scale_color_manual(values = clone_pal, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("CNV Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH1i-02-T1_CRv5_clone_pal", 4, 4)

### Assign cells to clones 
meta_IDH1i_02 = meta_IDH1i_02 %>% mutate(clone_manual = case_when(clone %in% c("gain_6_10", "gain_10") ~ "gain_6_10_clone",
                                                                  seurat_clusters %in% c(4, 6, 1, 9, 8, 7) ~ "NRAS_clone",
                                                                  seurat_clusters %in% c(0,2,3,5,10,11,12,13,14,15,16, 17) ~ "parent_clone"))


p_all_umap = ggplot(meta_IDH1i_02,  aes(x=IDH1i_02_UMAP_1, y=IDH1i_02_UMAP_2, color = clone_manual) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point( size=umap_pt_size) +  scale_color_manual(values = pal_clones_manual, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH1i-02-T1_clone_assignment", 4, 4)

NRAS_plot = meta_IDH1i_02 %>% dplyr::filter(NRAS_G12 %in% c("WT", "MT")) %>% 
    mutate(NRAS_G12 = factor(NRAS_G12, levels = c("WT", "MT"))) %>% 
    mutate(clone_manual = factor(clone_manual, levels = c("parent_clone", "NRAS_clone", "gain_6_10_clone")) )
p_stack = ggplot(NRAS_plot, aes(x = clone_manual, fill = NRAS_G12)) + 
    geom_bar(position = "fill", stat = "count") + scale_fill_manual(values = NRAS_pal, na.value = makeTransparent("grey"))+ theme_classic()  + 
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Clone") + ylab("Fraction Genotyped Cells") 
save_plot(p_stack, "stack_IDH1i-02-T1_NRAS_per_clone", 2.5, 4)

IDH_plot = meta_IDH1i_02 %>% dplyr::filter(IDH1_R132 %in% c("WT", "MT")) %>% 
    mutate(IDH1_R132 = factor(IDH1_R132, levels = c("WT", "MT"))) %>% 
    mutate(clone_manual = factor(clone_manual, levels = c("parent_clone", "NRAS_clone", "gain_6_10_clone")) )
p_stack = ggplot(IDH_plot, aes(x = clone_manual, fill = IDH1_R132)) + 
    geom_bar(position = "fill", stat = "count") + scale_fill_manual(values = IDH1_pal, na.value = makeTransparent("grey"))+ theme_classic()  + 
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Clone") + ylab("Fraction Genotyped Cells") 
save_plot(p_stack, "stack_IDH1i-02-T1_IDH1_per_clone", 2.5, 4)

table(meta_IDH1i_02$clone_manual, meta_IDH1i_02$NRAS_G12, useNA = "always")
table(meta_IDH1i_02$clone_manual, meta_IDH1i_02$IDH1_R132, useNA = "always")



## plot expression of surface markers 

genes = c("ANPEP", "ITGAM", "CD33")
for (module in genes){
  IDH1i_02_obj = calc_module_score_log10_scaled(IDH1i_02_obj, module, as.name(module))
}
meta_df_IDH1i_02_obj = IDH1i_02_obj[[]]

p_markers = list()
for (i in 1:length(genes)){
  meta_df_ordered_exp = meta_df_IDH1i_02_obj %>% arrange(.data[[paste0(genes[i], "_score")]])
p_markers[[i]] = ggplot(meta_df_ordered_exp, aes_string(x = "IDH1i_02_UMAP_1", y = "IDH1i_02_UMAP_2", color = paste0(genes[i], "_score")) )  + 
  geom_point(size = 0.01) + ggtitle(genes[i])  + theme_classic() + scale_colour_gradient(low = "gray80",
  high = muted("red")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + labs(color=' ') + xlab("") + ylab("")
        
}
leg = cowplot::get_legend(  p_markers[[1]] + theme(legend.position = "right"))
umap_markers = cowplot::plot_grid(plotlist = p_markers, nrow = 1, ncol = 3)
umap_markers_leg = cowplot::plot_grid(umap_markers, leg, nrow = 1, ncol = 2, rel_widths = c(3,0.5))
ggsave(plot = umap_markers_leg, filename = paste0(outdir_all, "figures_integrated_object/", "umap_scaled_exp_IDH1i_02_marker_genes_2_recolor.png", sep = ""),width = 6.5,height = 2) #,dpi = 500)
ggsave(plot = umap_markers_leg, filename = paste0(outdir_all, "figures_integrated_object/", "umap_scaled_exp_IDH1i_02_marker_genes_2_recolor.pdf", sep = ""),width = 6.5,height = 2) #,dpi = 500)


### MODULE GENE EXPRESSION BY CLONE 

###### boxplot of HSC and Imm Mye modules by clone 
lin_scores = c("vanGalen19_AML_HSC.like", "vanGalen19_AML_Myeloid.like", "S Phase Score") #, "SPI1.GFI1", "HLF.ZFP36L2", "vanGalen19_AML_Promono.like")
# score_vars = paste0(lin_scores, "_module_score_log2UMI")
score_vars = c("vanGalen19_AML_HSC.like_seurat.score1", "vanGalen19_AML_Myeloid.like_seurat.score1", "vanGalen19_AML_Myeloid.like_seurat.score1",  "S.Score")
y_labs_scores = c("HSC-like AML Score", "Myeloid-like AML Score", "GMP-like AML Score", "Cell Cycle S Phase Score")
my_comparisons <- list( c("parent_clone", "NRAS_clone") )
# stat_compare_means(comparisons = my_comparisons, size = 2.5) 

meta_IDH1i_02_subset = subset(meta_IDH1i_02, clone_manual %in% c("parent_clone", "NRAS_clone"))
meta_IDH1i_02_subset$clone_manual = factor(meta_IDH1i_02_subset$clone_manual, levels = c("parent_clone", "NRAS_clone"))
p_boxplot_lin_scores = list()
p_boxplot_lin_scores_no_leg = list()
  meta_IDH1i_02_subset[[paste0("vanGalen19_AML_HSC.like_seurat.score1", "_scaled")]] = scales::rescale(meta_IDH1i_02_subset[["vanGalen19_AML_HSC.like_seurat.score1"]], to = c(0,1))
  meta_IDH1i_02_subset[[paste0("vanGalen19_AML_Myeloid.like_seurat.score1", "_scaled")]] = scales::rescale(meta_IDH1i_02_subset[["vanGalen19_AML_Myeloid.like_seurat.score1"]], to = c(0,1))
  meta_IDH1i_02_subset[[paste0("S.Score", "_scaled")]] = scales::rescale(meta_IDH1i_02_subset[["S.Score"]], to = c(0,1))

for (i in 1:length(score_vars)){
p_boxplot_lin_scores[[i]] = ggplot(meta_IDH1i_02_subset,  aes_string(x="clone_manual", y= paste0(score_vars[i], "_scaled")) )+ geom_boxplot(aes(fill = clone_manual), outlier.shape = NA) + 
 scale_fill_manual(values = pal_clones_manual, limits = force)   + ylab(lin_scores[i]) + xlab("Clone") + ylab(y_labs_scores[i]) + 
 stat_compare_means(comparisons = my_comparisons, size = 2.5, label.y = 0.9) + 
 theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1, vjust = 1)) 
 p_boxplot_lin_scores_no_leg[[i]]  = p_boxplot_lin_scores[[i]] +  theme(legend.position = "none") + scale_x_discrete(labels=c("parent_clone" = "Parent Clone", "NRAS_clone" = "NRAS Clone"))
}

# p_boxplot_lin_scores_legend = get_legend(  p_boxplot_lin_scores_no_leg[[1]] )
cow_lin_bp = cowplot::plot_grid(plotlist = p_boxplot_lin_scores_no_leg, nrow = 1, ncol = 3)
save_plot_simple(cow_lin_bp, "IDH1i-02_pre_clone_lin_scores_bp", 5, 3)


### DIFFERENTIAL GENE EXPRESSION

# MPPs only 
NRAS_clone_MPP = meta_IDH1i_02 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("NRAS_clone"))) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH1i_02 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
NRAS_clone_v_parent_clone_MPP = FindMarkers(IDH1i_02_obj, ident.1 = rownames(NRAS_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

NRAS_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
NRAS_clone_v_parent_clone_MPP['gene'] = rownames(NRAS_clone_v_parent_clone_MPP)
write.csv(NRAS_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_NRAS_clone_v_parent_clone_MPP_logFC0.csv"), quote = FALSE)
NRAS_clone_v_parent_clone_MPP = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_NRAS_clone_v_parent_clone_MPP_logFC0.csv"), row.names = 1)

outdir = outdir_all

doFGSEA(NRAS_clone_v_parent_clone_MPP, "NRAS_clone_v_parent_clone_MPP", "NRAS Clone v Parent Clone MPPs")



saveRDS(IDH1i_02_obj, "/home/aifantis_collab/IDH1i-02_subset_clone_assigned.rds")
saveRDS(IDH2i_01_obj, "/home/aifantis_collab/IDH2i-01_subset_clone_assigned.rds")



meta_nras = IDH1i_02_obj[[]]
meta_nras_export = meta_nras %>% dplyr::select("IDH1i_02_UMAP_1", "IDH1i_02_UMAP_2", "ANPEP_score", "ITGAM_score", "CD33_score", "clone_manual", "HOXA9...", "HIF1A...", "RELA...", "RELB...", "STAT3...", "STAT1...", "SPI1...", "CEBPA...")


