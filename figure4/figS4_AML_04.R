## Figure 4 

source("../plotting_helpers.R")

## Figure params 
umap_pt_size = 0.2

IDH2_pal = c("#5c315f", "#c28c09", makeTransparent("grey"))
names(IDH2_pal) = c("MT", "WT", "NA")
SRSF2_pal = c("#0C4F4E", "#c28c09", makeTransparent("grey"))
names(SRSF2_pal) = c("MT", "WT", "NA")

## Read in object 

IDH2i_01_obj = readRDS("AML-04_object.rds")
meta_IDH2i_01 = IDH2i_01_obj[[]]
meta_IDH2i_01 = meta_IDH2i_01 %>% mutate(cell_type_simple = case_when(manual_annotation_v3_l1_final %in% c("Immature Mono") ~ "MPP", 
                                                                      TRUE ~ as.character(manual_annotation_v3_l1_final)))


## plot UMAPs colored by celltype and genotype 
p_all_umap = ggplot(meta_IDH2i_01,  aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point( aes(color = (cell_type_simple)), size=umap_pt_size)+ scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    theme_umap() + ggtitle("AML-05 UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i-01_T1", 4, 4)


p_all_umap = ggplot(meta_IDH2i_01, aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = SRSF2_P95) ) + 
    geom_point(data = subset(meta_IDH2i_01, (is.na(SRSF2_P95))),  size=umap_pt_size) +  scale_color_manual(values = SRSF2_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH2i_01, (SRSF2_P95 == "WT")), size=umap_pt_size) +  scale_color_manual(values = SRSF2_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH2i_01, (SRSF2_P95 == "MT")), size=umap_pt_size) +  scale_color_manual(values = SRSF2_pal, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("SRSF2 P95")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i-01-T1_CRv5_color_SRSF2_P95", 4, 4)

p_all_umap = ggplot(meta_IDH2i_01, aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = IDH2_R172) ) + 
    geom_point(data = subset(meta_IDH2i_01, (is.na(IDH2_R172))),  size=umap_pt_size) +  scale_color_manual(values = IDH2_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH2i_01, (IDH2_R172 == "WT")), size=umap_pt_size) +  scale_color_manual(values = IDH2_pal, na.value = makeTransparent("grey"))  +
    geom_point(data = subset(meta_IDH2i_01, (IDH2_R172 == "MT")), size=umap_pt_size) +  scale_color_manual(values = IDH2_pal, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("IDH2 R172")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i-01-T1_CRv5_color_IDH2_R172", 4, 4)

cnv_pal_2 = c("steelblue3", "yellow", makeTransparent("grey"))
 names(cnv_pal_2) = c("gain_14", "gain_8", "normal")

p_all_umap = ggplot(meta_IDH2i_01, aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = clone) ) + 
    geom_point(  size=umap_pt_size) +  scale_color_manual(values = cnv_pal_2, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("CNV Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i-01-T1_CRv5_cnv_clone", 4, 4)



# assign cells to clones
meta_IDH2i_01 = meta_IDH2i_01 %>% mutate(clone_manual = case_when(seurat_clusters %in% c(10,4) ~ "parent_clone",
                                                                  seurat_clusters %in% c(0,1,2,3,4,5,7,9,10,11,12) ~ "SRSF2_subclone",
                                                                  seurat_clusters %in% c(6) ~ "mono_1",
                                                                  seurat_clusters %in% c(8) ~ "mono_2"))

clone_pal_IDH2 = c("#5c315f", "#0C4F4E")
names(clone_pal_IDH2) = c("parent_clone", "SRSF2_subclone")

## stacked barplot of mutant cells per clone 

p_all_umap = ggplot(meta_IDH2i_01,  aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = clone_manual) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point( size=umap_pt_size) +  scale_color_manual(values = clone_pal_IDH2, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i_01-T1_clone_assignment", 4, 4)

SRSF2_plot = meta_IDH2i_01 %>% dplyr::filter(SRSF2_P95 %in% c("WT", "MT")) %>% 
dplyr::filter(!is.na(clone_manual)) %>% 
    mutate(SRSF2_P95 = factor(SRSF2_P95, levels = c("WT", "MT"))) %>% 
    mutate(clone_manual = case_when(clone_manual %in% c("mono_1", "mono_2") ~ "mono", TRUE ~ clone_manual)) %>% 
    mutate(clone_manual = factor(clone_manual, levels = c("parent_clone", "SRSF2_subclone", "mono")) )
p_stack = ggplot(SRSF2_plot, aes(x = clone_manual, fill = SRSF2_P95)) + 
    geom_bar(position = "fill", stat = "count") + scale_fill_manual(values = SRSF2_pal, na.value = makeTransparent("grey"))+ theme_classic()  + 
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Clone") + ylab("Fraction Genotyped Cells") 
save_plot(p_stack, "stack_IDH2i_01-T1_SRSF2_per_clone_nomono", 2.5, 4)

IDH2_plot = meta_IDH2i_01 %>% dplyr::filter(IDH2_R172 %in% c("WT", "MT")) %>% 
    mutate(IDH2_R172 = factor(IDH2_R172, levels = c("WT", "MT"))) %>% 
    mutate(clone_manual = case_when(clone_manual %in% c("mono_1", "mono_2") ~ "mono", TRUE ~ clone_manual)) %>% 
    mutate(clone_manual = factor(clone_manual, levels = c("parent_clone", "SRSF2_subclone", "mono")) )
p_stack = ggplot(IDH2_plot, aes(x = clone_manual, fill = IDH2_R172)) + 
    geom_bar(position = "fill", stat = "count") + scale_fill_manual(values = IDH2_pal, na.value = makeTransparent("grey"))+ theme_classic()  + 
    theme( axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Clone") + ylab("Fraction Genotyped Cells") 
save_plot(p_stack, "stack_IDH2i_01-T1_IDH2_per_clone_nomono", 2.5, 4)

table(meta_IDH2i_01$clone_manual, meta_IDH2i_01$SRSF2_P95, useNA = "always")
table(meta_IDH2i_01$clone_manual, meta_IDH2i_01$IDH2_R172, useNA = "always")

# color mono 1 mono 2 
clone_pal_IDH2_mono = c("#5c315f", "#0C4F4E", "#a77aa0", "#7aa0a7")
names(clone_pal_IDH2_mono) = c("parent_clone", "SRSF2_subclone", "mono_1", "mono_2")

p_all_umap = ggplot(meta_IDH2i_01,  aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = clone_manual) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point( size=umap_pt_size) +  scale_color_manual(values = clone_pal_IDH2_mono, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i_01-T1_clone_assignment_mono", 4, 4)


clone_pal_IDH2_mono = c("#5c315f", "#0C4F4E", "#70A9A1")
names(clone_pal_IDH2_mono) = c("parent_clone", "SRSF2_subclone", "mono")

p_all_umap = ggplot(meta_IDH2i_01,  aes(x=IDH2i_01_UMAP_1, y=IDH2i_01_UMAP_2, color = clone_manual) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point( size=umap_pt_size) +  scale_color_manual(values = clone_pal_IDH2_mono, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    ggtitle("Clones")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_IDH2i_01-T1_clone_assignment_nomono", 4, 4)





###### boxplot of HSC and Imm Mye modules by clone 
library(ggpubr)
lin_scores = c("vanGalen19_AML_HSC.like", "vanGalen19_AML_Myeloid.like", "S Phase Score") #, "SPI1.GFI1", "HLF.ZFP36L2", "vanGalen19_AML_Promono.like")
# score_vars = paste0(lin_scores, "_module_score_log2UMI")
score_vars = c("vanGalen19_AML_HSC.like_seurat.score1", "vanGalen19_AML_Myeloid.like_seurat.score1", "TAL1.HSF1_seurat.score1", "vanGalen19_AML_Progenitor.like_seurat.score1", "vanGalen19_AML_Promono.like_seurat.score1", "S.Score", "G2M.Score")
y_labs_scores = c("HSC-like AML Score", "Myeloid-like AML Score", "Ery-like Score", "AML Progenitor-like Score", "AML Promono-like", "Cell Cycle S Phase Score", "Cell Cycle G2M Phase Score")
my_comparisons <- list( c("parent_clone", "SRSF2_subclone") )
# stat_compare_means(comparisons = my_comparisons, size = 2.5) 

meta_IDH2i_01_subset = subset(meta_IDH2i_01, clone_manual %in% c("parent_clone", "SRSF2_subclone"))
meta_IDH2i_01_subset$clone_manual = factor(meta_IDH2i_01_subset$clone_manual, levels = c("parent_clone", "SRSF2_subclone"))
p_boxplot_lin_scores = list()
p_boxplot_lin_scores_no_leg = list()
for (i in 1:length(score_vars)){
p_boxplot_lin_scores[[i]] = ggplot(meta_IDH2i_01_subset,  aes_string(x="clone_manual", y=score_vars[i]) ) + geom_boxplot(aes(fill = clone_manual), outlier.shape = NA) + 
 scale_fill_manual(values = clone_pal_IDH2, limits = force)   + ylab(lin_scores[i]) + xlab("Clone") + ylab(y_labs_scores[i]) + 
 stat_compare_means(comparisons = my_comparisons, size = 2.5) + 
 theme_classic() + theme(axis.text.x=element_text(angle=45,hjust=1, vjust = 1)) 
 p_boxplot_lin_scores_no_leg[[i]]  = p_boxplot_lin_scores[[i]] +  theme(legend.position = "none") + scale_x_discrete(labels=c("parent_clone" = "Parent Clone", "SRSF2_subclone" = "SRSF2 Clone"))
}

cow_lin_bp = cowplot::plot_grid(plotlist = p_boxplot_lin_scores_no_leg, nrow = 1, ncol = length(score_vars))
ggsave(plot = cow_lin_bp, width = 8, height = 3, filename = paste0(outdir_all, "figures_integrated_object/IDH2i-01_pre_clone_lin_scores_bp.png"))
ggsave(plot = cow_lin_bp, width = 8, height = 3, filename = paste0(outdir_all, "figures_integrated_object/IDH2i-01_pre_clone_lin_scores_bp.pdf"))

### differential expression between IDH2 and SRSF2 clone
SRSF2_clone = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% tibble::column_to_rownames('cb')
parent_clone = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone), ident.2 = rownames(parent_clone), logfc.threshold = 0)

#SRSF2_clone_v_parent_clone %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone['gene'] = rownames(SRSF2_clone_v_parent_clone)
write.csv(SRSF2_clone_v_parent_clone, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_logFC0.csv"), quote = FALSE)
SRSF2_clone_v_parent_clone = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_logFC0.csv"), row.names = 1)

# MPPs only 
SRSF2_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% dplyr::filter(clone %in% c("normal", "gain_8")) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone_MPP = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

#SRSF2_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone_MPP['gene'] = rownames(SRSF2_clone_v_parent_clone_MPP)
write.csv(SRSF2_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_excludegain8_14.csv"), quote = FALSE)
SRSF2_clone_v_parent_clone_MPP = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_excludegain8_14.csv"), row.names = 1)

##### try different combination of excluding/including CNVs 

# MPPs only 
SRSF2_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% dplyr::filter(clone %in% c("normal")) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone_MPP = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

#SRSF2_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone_MPP['gene'] = rownames(SRSF2_clone_v_parent_clone_MPP)
write.csv(SRSF2_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2excludegain8_14_allidh.csv"), quote = FALSE)
SRSF2_clone_v_parent_clone_MPP = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2excludegain8_14_allidh.csv"), row.names = 1)

fgsea_res_SRSF2 = doFGSEA(SRSF2_clone_v_parent_clone_MPP, "srsf2excludegain8_14_allidh", "SRSF2 Clone v Parent Clone MPPs")

fgsea_res_SRSF2_plot = fgsea_res_SRSF2 %>% dplyr::filter(pathway %in% c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA", "HALLMARK_DNA_REPAIR")) %>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% mutate(pathway = gsub("HALLMARK", "", pathway))
p_gsea <- ggplot(fgsea_res_SRSF2_plot, aes(x = NES, y = fct_reorder(pathway, NES))) + 
               geom_point(aes(size = abs(NES), color = padj)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        xlab("Normalized Enrichment Score") + 
        ylab("Hallmark Pathway") + # facet_grid(.~type, scales="free") +
        ggtitle("Pathways Enriched in SRSF2 Subclone vs Parent Clone")
save_plot_simple(p_gsea, "gsea_dotplot_SRSF2_subset_srsf2excludegain8_14_allidh", 7, 4)


# MPPs only 
SRSF2_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% dplyr::filter(clone %in% c("normal", "gain_8")) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone_MPP = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

SRSF2_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone_MPP['gene'] = rownames(SRSF2_clone_v_parent_clone_MPP)
write.csv(SRSF2_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2exclude14_allidh.csv"), quote = FALSE)

fgsea_res_SRSF2 = doFGSEA(SRSF2_clone_v_parent_clone_MPP, "srsf2exclude14_allidh", "SRSF2 Clone v Parent Clone MPPs")

fgsea_res_SRSF2_plot = fgsea_res_SRSF2 %>% dplyr::filter(pathway %in% c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA", "HALLMARK_DNA_REPAIR")) %>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% mutate(pathway = gsub("HALLMARK", "", pathway))
p_gsea <- ggplot(fgsea_res_SRSF2_plot, aes(x = NES, y = fct_reorder(pathway, NES))) + 
               geom_point(aes(size = abs(NES), color = padj)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        xlab("Normalized Enrichment Score") + 
        ylab("Hallmark Pathway") + # facet_grid(.~type, scales="free") +
        ggtitle("Pathways Enriched in SRSF2 Subclone vs Parent Clone")
save_plot_simple(p_gsea, "gsea_dotplot_srsf2exclude14_allidh", 7, 4)



# MPPs only 
SRSF2_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% dplyr::filter(clone %in% c("normal")) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(clone %in% c("normal"))%>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone_MPP = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

#SRSF2_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone_MPP['gene'] = rownames(SRSF2_clone_v_parent_clone_MPP)
write.csv(SRSF2_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2normalonly_idhnormalonly.csv"), quote = FALSE)


fgsea_res_SRSF2 = doFGSEA(SRSF2_clone_v_parent_clone_MPP, "srsf2normalonly_idhnormalonly", "SRSF2 Clone v Parent Clone MPPs")

fgsea_res_SRSF2_plot = fgsea_res_SRSF2 %>% dplyr::filter(pathway %in% c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA", "HALLMARK_DNA_REPAIR")) %>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% mutate(pathway = gsub("HALLMARK", "", pathway))
p_gsea <- ggplot(fgsea_res_SRSF2_plot, aes(x = NES, y = fct_reorder(pathway, NES))) + 
               geom_point(aes(size = abs(NES), color = padj)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        xlab("Normalized Enrichment Score") + 
        ylab("Hallmark Pathway") + # facet_grid(.~type, scales="free") +
        ggtitle("Pathways Enriched in SRSF2 Subclone vs Parent Clone")
save_plot_simple(p_gsea, "gsea_dotplot_srsf2normalonly_idhnormalonly", 7, 4)

# MPPs only 
SRSF2_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("SRSF2_subclone"))) %>% dplyr::filter(clone %in% c("gain_8")) %>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>% tibble::column_to_rownames('cb')
parent_clone_MPP = meta_IDH2i_01 %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( clone_manual%in% c("parent_clone"))) %>% dplyr::filter(clone %in% c("gain_8"))%>% dplyr::filter(manual_annotation_v3_l2_final %in% c("MPP")) %>%  tibble::column_to_rownames('cb')
SRSF2_clone_v_parent_clone_MPP = FindMarkers(IDH2i_01_obj, ident.1 = rownames(SRSF2_clone_MPP), ident.2 = rownames(parent_clone_MPP), logfc.threshold = 0)

SRSF2_clone_v_parent_clone_MPP %>% arrange(avg_log2FC)
SRSF2_clone_v_parent_clone_MPP['gene'] = rownames(SRSF2_clone_v_parent_clone_MPP)
write.csv(SRSF2_clone_v_parent_clone_MPP, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2only8_idhonly8.csv"), quote = FALSE)


fgsea_res_SRSF2 = doFGSEA(SRSF2_clone_v_parent_clone_MPP, "srsf2only8_idhonly8", "SRSF2 Clone v Parent Clone MPPs")

fgsea_res_SRSF2_plot = fgsea_res_SRSF2 %>% dplyr::filter(pathway %in% c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA", "HALLMARK_DNA_REPAIR")) %>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% mutate(pathway = gsub("HALLMARK", "", pathway))
p_gsea <- ggplot(fgsea_res_SRSF2_plot, aes(x = NES, y = fct_reorder(pathway, NES))) + 
               geom_point(aes(size = abs(NES), color = padj)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        xlab("Normalized Enrichment Score") + 
        ylab("Hallmark Pathway") + # facet_grid(.~type, scales="free") +
        ggtitle("Pathways Enriched in SRSF2 Subclone vs Parent Clone")
save_plot_simple(p_gsea, "gsea_dotplot_srsf2only8_idhonly8", 7, 4)




outdir = outdir_all

doFGSEA(SRSF2_clone_v_parent_clone_MPP, "SRSF2_clone_v_parent_clone_MPP", "SRSF2 Clone v Parent Clone MPPs")


## read in post treatment objects
IDH2i_01_obj_T2 = readRDS("all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_081722_subset_IDH2i_01_obj_T2_updated_numbat.rds")
meta_IDH2i_01_T2 = IDH2i_01_obj_T2[[]]
IDH2i_01_obj_T3 = readRDS("all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_081722_subset_IDH2i_01_obj_T3_numbat_updated.rds")
meta_IDH2i_01_T3 = IDH2i_01_obj_T3[[]]


mut_stack_df = bind_rows(meta_IDH2i_01, meta_IDH2i_01_T2, meta_IDH2i_01_T3)
mut_stack_df$cell_type = mut_stack_df$manual_annotation_v3_l1_final
mut_stack_df = mut_stack_df %>% dplyr::filter(cell_type %in% c("HSC", "MPP", "MEP", "Ery", "CD14+ Mono" , "CD16+ Mono", "cDC", "Immature Mono") ) %>% 
                                                                          mutate(cell_type_broad = case_when(cell_type %in% c("HSC", "MPP", "MEP", "Ery") ~ "HSPC",
                                                                          cell_type %in% c("CD14+ Mono" , "CD16+ Mono", "cDC", "Immature Mono") ~ "Mono/DC"
                                                                         
                                                                          )) %>% mutate(cell_type_broad = factor(cell_type_broad, levels = c( "HSPC", "Mono/DC" ))) 
mut_stack_df_SRSF2 = subset(mut_stack_df, SRSF2_P95 %in% c("MT", "WT"))
p_stacked_mt = ggplot(mut_stack_df_SRSF2, aes(fill=factor(SRSF2_P95, levels = c("WT", "MT")),  x=Timepoint)) + 
    geom_bar(position="fill", stat="count") + scale_fill_manual(values = SRSF2_pal, na.value = makeTransparent("grey"))  +
    ylab("Fraction") + xlab('')+ theme_bw() + facet_grid(. ~ cell_type_broad, scales = "free", space='free')  + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
theme(legend.position='bottom') + labs(fill=' ') + scale_x_discrete(breaks=c("T1","T2","T3"),
        labels=c("Diagnosis", "Persistent AML (D28)", "Persistent AML (D84)")) + 
ggtitle("")
# ct_legend = cowplot::get_legend(p_stacked_ct_normal + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 2, title = "Cell Type")))
ggsave(plot = p_stacked_mt, width = 2, height = 3, filename = paste0(outdir_all, "figures_integrated_object/", "stacked_barplot_IDH2i-01_celltype_SRSF2.png"))
ggsave(plot = p_stacked_mt, width = 2, height = 3, filename = paste0(outdir_all, "figures_integrated_object/", "stacked_barplot_IDH2i-01_celltype_SRSF2.pdf"))

mut_stack_df_IDH2 = subset(mut_stack_df, IDH2_R172 %in% c("MT", "WT"))
p_stacked_mt = ggplot(mut_stack_df_IDH2, aes(fill=factor(IDH2_R172, levels = c("WT", "MT")),  x=Timepoint)) + 
    geom_bar(position="fill", stat="count") + scale_fill_manual(values = IDH2_pal, na.value = makeTransparent("grey"))  +
    ylab("Fraction") + xlab('')+ theme_bw() + facet_grid(. ~ cell_type_broad, scales = "free", space='free')  + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
theme(legend.position='bottom') + labs(fill=' ') + scale_x_discrete(breaks=c("T1","T2","T3"),
        labels=c("Diagnosis", "Persistent AML (D28)", "Persistent AML (D84)")) + 
ggtitle("")
# ct_legend = cowplot::get_legend(p_stacked_ct_normal + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 2, title = "Cell Type")))
ggsave(plot = p_stacked_mt, width = 2, height = 3, filename = paste0(outdir_all, "figures_integrated_object/", "stacked_barplot_IDH2i-01_celltype_IDH2.png"))
ggsave(plot = p_stacked_mt, width = 2, height = 3, filename = paste0(outdir_all, "figures_integrated_object/", "stacked_barplot_IDH2i-01_celltype_IDH2.pdf"))





########
## Blast co-mutation groups -- co-mutation modules 

### assign all dx blasts to a co-mutation group 
source("/home/seurat/source_files/heatmap_DEGS.R")

meta_all = int_obj_all[[]]


### update donor.id 
names_dict = list("IDH1i-01" = "AML-01", "IDH1i-02" = "AML-02", "IDH1i-03" = "AML-03",
                  "IDH2i-01" = "AML-04", "IDH2i-02" = "AML-05", "IDH2i-03" = "AML-06",
                  "AML0048" = "AML-07", "AML2123" = "AML-08", "AML1529" = "AML-09",
                  "AML4897" = "AML-10", "AML0875" = "AML-11", "AML1371" = "AML-12",
                  "AML4340" = "AML-13", "AML0024" = "AML-14", "AML1896" = "AML-15", 
                  "CH-01" = "CH-01", "CH-05" = "CH-02", 
                  "healthy_MSK_01" = "Donor-01-PB", "healthy_MSK_02" = "Donor-02-PB", 
                  "healthy_MSK_03" = "Donor-03-BM", "healthy_MSK_04" = "Donor-04-BM", 
                  "BM_healthy_NYU_1" = "Donor-05-BM", "BM_healthy_NYU_2" = "Donor-06-BM",
                  "BM_healthy_NYU_3" = "Donor-07-BM", "BM_healthy_NYU_4" = "Donor-08-BM",
                  "BM_healthy_NYU_5" = "Donor-09-BM", "BM_healthy_NYU_6" = "Donor-10-BM",
                  "BM_healthy_NYU_7" = "Donor-11-BM", "BM_healthy_NYU_8" = "Donor-12-BM"
                  )


meta_all = meta_all %>% mutate(donor.id.original = donor.id) 
for (i in 1:length(names_dict)){meta_all$donor.id <- replace(meta_all$donor.id, meta_all$donor.id == names(names_dict[i]), as.character(names_dict[i]))}

### rename immature mono as MPP 
ct_dict = list("Immature Mono" = "MPP")
meta_all = meta_all %>% mutate(manual_annotation_v3_l1_final.original = manual_annotation_v3_l1_final) 
for (i in 1:length(ct_dict)){meta_all$manual_annotation_v3_l1_final <- replace(meta_all$manual_annotation_v3_l1_final, meta_all$manual_annotation_v3_l1_final == names(ct_dict[i]), as.character(ct_dict[i]))}
meta_all$manual_annotation_v3_l1_final = droplevels(meta_all$manual_annotation_v3_l1_final)

meta_all0 = meta_all %>% dplyr::filter(!(donor.id == "CH-06"))

meta_all_dx = meta_all0  %>% filter(institute_disease %in% c("MSK_AML"), stage_simple == "diagnosis")



# meta_IDH1i_02 %>% select(clone_manual) %>% write.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/IDH1i_02_clone_labels.csv"), quote = FALSE)
# meta_IDH2i_01 %>% select(clone_manual) %>%  write.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/IDH2i_01_clone_labels_test.csv"), quote = FALSE)

NRAS_clone = meta_IDH1i_02 %>% filter(clone_manual == "NRAS_clone") 
IDH1_clone = meta_IDH1i_02 %>% filter(clone_manual == "parent_clone") 

SRSF2_clone = meta_IDH2i_01 %>% filter(clone_manual == "SRSF2_subclone") %>% filter(clone == "gain_8") 
IDH2_clone = meta_IDH2i_01 %>% filter(clone_manual == "parent_clone") %>% filter(clone == "gain_8") 



meta_all_dx = meta_all_dx %>% tibble::rownames_to_column("cb") %>% mutate(co_mut_group = case_when(cb %in% rownames(NRAS_clone) ~ "NRAS",
                                                              cb %in% rownames(IDH1_clone) ~ "IDH1",
                                                              cb %in% rownames(SRSF2_clone) ~ "SRSF2",
                                                              cb %in% rownames(IDH2_clone) ~ "IDH2",
                                                              donor.id.original %in% c("IDH1i-01", "IDH1i-03") ~ "NPM1", 
                                                              donor.id.original %in% c("IDH2i-02", "IDH2i-03") ~ "SRSF2", 
                                                              TRUE ~ "other"
)) %>% tibble::column_to_rownames("cb")  %>% filter(manual_annotation_v3_l1_final %in% c("HSC","MPP","GMP", "MEP", "Ery" ))%>% filter(!(co_mut_group %in% c("other" )))



# co_mut_subset_obj = subset(int_obj_all, cells = rownames(meta_all_dx))



meta_all_dx_MPP = meta_all_dx %>% dplyr::filter(manual_annotation_v3_l1_final %in% "MPP")

co_mut_subset_obj_MPP = subset(int_obj_all, cells = rownames(meta_all_dx_MPP))

co_mut_subset_obj_MPP = AddMetaData(co_mut_subset_obj_MPP, meta_all_dx_MPP['co_mut_group'])

Idents(co_mut_subset_obj_MPP) = "co_mut_group"
DefaultAssay(co_mut_subset_obj_MPP) = "RNA"


## downsample blasts by patient 
co_mut_subset_obj_MPP$donor_group = paste0(co_mut_subset_obj_MPP$donor.id, "_", co_mut_subset_obj_MPP$co_mut_group)
Idents(co_mut_subset_obj_MPP) = "donor_group"
MPP.small <- subset(co_mut_subset_obj_MPP, downsample = 300)


Idents(MPP.small) = "co_mut_group"
MPP_small_markers_MPP = FindAllMarkers(MPP.small , logfc.threshold = 0.25)

MPP.small <- ScaleData(object = MPP.small, features = rownames(MPP.small))

MPP_small_markers_MPP_0 = MPP_small_markers_MPP %>% dplyr::filter(!(stringr::str_detect(gene, "^RPL")))%>% dplyr::filter( !(stringr::str_detect(gene, "^RPS")))%>% dplyr::filter( !(stringr::str_detect(gene, "^MT-")))
top10<-MPP_small_markers_MPP_0 %>% group_by(cluster) %>% top_n(30, avg_log2FC)
MPP_small_markers_MPP_0 %>%   write.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/co_mut_modules_MSK_cohort_MPP_logFC0.25_figS4_idh2srsf2gain8only.csv"), quote = FALSE)


co_mut_subset_obj_MPP <- ScaleData(object = co_mut_subset_obj_MPP, features = rownames(co_mut_subset_obj_MPP))

co_mut_group_markers_MPP_0 = co_mut_group_markers_MPP %>% dplyr::filter(!(stringr::str_detect(gene, "^RPL")))%>% dplyr::filter( !(stringr::str_detect(gene, "^RPS")))%>% dplyr::filter( !(stringr::str_detect(gene, "^MT-")))

co_mut_group_markers_MPP_0 %>%   write.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/co_mut_modules_MSK_cohort_MPP_logFC0.25_figS4_idh2srsf2gain8only.csv"), quote = FALSE)
co_mut_group_markers_MPP_0 = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/co_mut_modules_MSK_cohort_MPP_logFC0.25_figS4_idh2srsf2gain8only.csv"), row.names = 1)

png(width=900, height=1500, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_MSK_MPP.small_figS4_idh2srsf2gain8only.png") )   
p = plot_heatmap_MS(dataset = MPP.small, 
              markers = gsub(" ", "", unique(as.character(top10$gene))),
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id","S.Score","G2M.Score"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660"), # color vector
                                 c("blue","white","red"),                            # Three-color gradient
                                 c("blue","white","red"))
)
dev.off()

TFs = scan("/home/seurat/patient_centric/TF_names_v_1.01.txt", what = character() ) 

MPP_small_markers_MPP_TFs = MPP_small_markers_MPP %>% dplyr::filter(gene %in% TFs) %>% dplyr::filter(!(stringr::str_detect(gene, "^RPL")))%>% dplyr::filter( !(stringr::str_detect(gene, "^RPS")))%>% dplyr::filter( !(stringr::str_detect(gene, "^MT-"))) %>% dplyr::filter( !(stringr::str_detect(gene, "^MALAT1"))) 
top10_TFs<-MPP_small_markers_MPP_TFs %>% group_by(cluster) %>% top_n(30, avg_log2FC)

png(width=900, height=1500, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_TFs_MSK_MPP.small_figS4_idh2srsf2gain8only.png") )
p = plot_heatmap_MS(dataset = MPP.small, 
              markers = gsub(" ", "", unique(as.character(top10_TFs$gene))),
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id","S.Score","G2M.Score"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660"), # color vector
                                 c("blue","white","red"),                            # Three-color gradient
                                 c("blue","white","red"))
)
dev.off()



subset_genes = c("NFE2", "HMGN3", "ITM2A", "CD74", "HLA-DRA", "HLA-DPB1", "HLA-DRB1", "LSP1", "CYTL1", "STMN1",  "ACTB", "RAC2", "ANGPT1", "CD69", # "PFN1", 
"SOX4", "SAT1", "MPO", "AZU1", "CFD", "ELANE", "CST3", "CSTA", "HOXB-AS3", "HCST", "SRGN", "S100A10", "VIM", "HOXA9", "HOXA10", "HOXA7",
"STAT3", "KDM5B", "CEBPA", 
"FTH1", "NR4A2", "NR4A1", "DDIT4",
"HMGA1", "OLIG1" )

MPP.small$donor.id.original = MPP.small$donor.id 
for (i in 1:length(names_dict)){MPP.small$donor.id  <- replace(MPP.small$donor.id , MPP.small$donor.id  == names(names_dict[i]), as.character(names_dict[i]))}


png(width=600, height=600, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_subset_MSK_figS4_idh2srsf2gain8only.png") )
p = plot_heatmap_MS(dataset = MPP.small,
              markers = unique(as.character(subset_genes)), ## can change to subset_genes for more restricted list 
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660") # color vector
                                )
)
dev.off()


png(width=600, height=1200, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_subset_labeled_MSK_figS4_idh2srsf2gain8only.png") )
p = plot_heatmap_MS_mark_annotation(dataset = MPP.small,
              markers = gsub(" ", "", unique(as.character(top10$gene))),
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660") # color vector
                                )
)
dev.off()


# as PDF 

pdf(width=6, height=6, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_subset_MSK_figS4_idh2srsf2gain8only.pdf") )
p = plot_heatmap_MS(dataset = MPP.small,
              markers = unique(as.character(subset_genes)), ## can change to subset_genes for more restricted list 
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660") # color vector
                                )
)
dev.off()
dev.off()


pdf(width=6, height=12, file= paste0(outdir_all, "figures_integrated_object/DGE_allmuts_blasts_scillus_heatmap_subset_labeled_MSK_figS4_idh2srsf2gain8only.pdf") )
p = plot_heatmap_MS_mark_annotation(dataset = MPP.small,
              markers = gsub(" ", "", unique(as.character(top10$gene))),
              sort_var = c("co_mut_group","donor.id"),
              anno_var = c("co_mut_group","donor.id"),
              anno_colors = list(co_mut_colors,                                             # RColorBrewer palette
                                 c("#191970", "#e17c7c", "#559cca", "#990000", "#1D7874", "#b46c4c","#FFBE7D", "#D37295",   "#8CD17D", 
"#B6992D", "#F1CE63",  "#86BCB6", "#E15759" , "#D4A6C8" ,"#9D7660") # color vector
                                )
)
dev.off()
dev.off()






meta_IDH2i_01_export = meta_IDH2i_01 %>% dplyr::select("IDH2i_01_UMAP_1", "IDH2i_01_UMAP_2", "clone_manual")
