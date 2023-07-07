## Figure 6 

source("../plotting_helpers.R")


CH_normal_obj_mono = readRDS("/home/aifantis_collab/data/all_combined_filtered_seurat_rpca_mt_rb_regressed_k10_081722_CH_normal_obj_mono.rds")
ch_got = meta_all0 %>% filter(donor.id %in% c("healthy_MSK_02", "healthy_MSK_01", "CH-01","CH-05"), manual_annotation_v3_l1_final %in% c("CD14+ Mono","CD16+ Mono", "Immature Mono", "cDC")) %>% select(IDH2_R140, SRSF2_P95)
CH_normal_obj_mono = AddMetaData(CH_normal_obj_mono, metadata=ch_got)
CH_normal_obj_mono_df = CH_normal_obj_mono[[]]


# meta_de = meta_query
umap_pt_size = 0.01
leg_size = 4
p_all_umap = ggplot(CH_normal_obj_mono_df,  aes(x=CH_UMAP_1, y=CH_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (institute_disease %in% c( "MSK_healthy"))), colour = "grey" , size=umap_pt_size) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (institute_disease %in% c( "MSK_CH"))), aes(color = as.factor(manual_annotation_v3_l1_final)) ,size=umap_pt_size) + #pch=21,
  scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
  theme_umap() + 
  ggtitle("CH UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_CH_normals", 4.3, 4.3)


p_all_umap = ggplot(CH_normal_obj_mono_df,  aes(x=CH_UMAP_1, y=CH_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point(data = subset(CH_normal_obj_mono_df, (institute_disease %in% c( "MSK_healthy"))), colour = "grey" , size=umap_pt_size) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (institute_disease %in% c( "MSK_CH"))), aes(color = as.factor(donor.id)) , #pch=21,
  size=umap_pt_size)+ scale_color_manual(values = c("#841c26", "#2C5784"), na.value = makeTransparent("grey"))  +
  theme_umap() + 
  ggtitle("CH UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_CH_patient", 4.3, 4.3)



ch_pt_pal = c("#841c26", "#2C5784")
names(ch_pt_pal) = c("CH-01", "CH-05")

p_all_umap = ggplot(CH_normal_obj_mono_df,  aes(x=CH_UMAP_1, y=CH_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (IDH2_R140 %in% c( "WT", NA))), colour = makeTransparent("grey") , size=umap_pt_size) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (IDH2_R140 %in% c( "MT"))), aes(color = as.factor(donor.id)) , #pch=21, color = "white", 
  size=umap_pt_size) + scale_color_manual(values = ch_pt_pal, na.value = makeTransparent("grey"))  +
  theme_umap() + 
  ggtitle("IDH2 R140 MT Cells ")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_CH_IDH2_R140", 4.3, 4.3)



p_all_umap = ggplot(CH_normal_obj_mono_df,  aes(x=CH_UMAP_1, y=CH_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (SRSF2_P95 %in% c( "WT", NA))), colour = makeTransparent("grey") , size=umap_pt_size) + 
  geom_point(data = subset(CH_normal_obj_mono_df, (SRSF2_P95 %in% c( "MT"))), aes(color = as.factor(donor.id)) , #pch=21, color = "white", 
  size=umap_pt_size)+ scale_color_manual(values = ch_pt_pal, na.value = makeTransparent("grey"))  +
  theme_umap() + 
  ggtitle("SRSF2 P95 MT Cells ")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = leg_size)))
save_plot(p_all_umap, "umap_CH_SRSF2_P95", 4.3, 4.3)



ch_mut_df = subset(CH_normal_obj_mono_df, donor.id %in% c("CH-01", "CH-05", "CH-02"))
ch_mut_df = ch_mut_df %>% mutate(ch_mut = case_when(IDH2_R140 == "MT" & SRSF2_P95 == "MT" ~ "IDH2 SRSF2 MUT",
                                                    IDH2_R140 == "MT" ~ "IDH2 MUT", 
                                                    SRSF2_P95 == "MT" ~ "SRSF2 MUT", 
                                                    IDH2_R140 == "WT" & SRSF2_P95 == "WT" ~ "IDH2 SRSF2 WT",
                                                    IDH2_R140 == "WT" ~ "IDH2 WT", 
                                                    SRSF2_P95 == "WT" ~ "SRSF2 WT", 
                                                    TRUE ~ "NA") )  %>% mutate(ch_mut = factor(ch_mut, levels = c("IDH2 SRSF2 MUT" , "IDH2 MUT", "SRSF2 MUT", "IDH2 SRSF2 WT", "IDH2 WT","SRSF2 WT", "NA" )))

# ch_mut_pal = c("#c89700","#5c315f" ,"#0c4f4e" , "#c1ccca", "#d3cbc7", "#aebcc4")
# names(ch_mut_pal) = c("IDH2 SRSF2 MUT", "IDH2 MUT", "SRSF2 MUT", "SRSF2 WT", "IDH2 WT","IDH2 SRSF2 WT")      

ch_mut_pal = c("#0c4f4e","#5c315f" ,"#AA510E" , "#c1ccca", "#d3cbc7", "#aebcc4") # EFA48B
names(ch_mut_pal) = c("IDH2 SRSF2 MUT", "IDH2 MUT", "SRSF2 MUT", "SRSF2 WT", "IDH2 WT","IDH2 SRSF2 WT")      


df_ch_mut_plot = subset(ch_mut_df, !(ch_mut %in% c("NA"))) 
df_ch_mut_plot$manual_annotation_v3_l1_final = droplevels(df_ch_mut_plot$manual_annotation_v3_l1_final)
p_stacked_ct = df_ch_mut_plot %>% group_by(donor.id, manual_annotation_v3_l1_final) %>% count(ch_mut) %>% 
  mutate(prop = n / sum(n)) %>%
  mutate(column_total = sum(n)) %>% ungroup() %>% 
ggplot( aes(fill=ch_mut,  x=manual_annotation_v3_l1_final, y = prop)) + 
    geom_col(position = position_fill(reverse = TRUE)
    # stat="count"
    ) + scale_fill_manual(values = ch_mut_pal, na.value = makeTransparent("grey"))  +
    ylab("Fraction of Genotyped Cells") + xlab('Cell Type')+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='bottom')  + facet_wrap(.~donor.id) + coord_flip() + 
guides(fill=guide_legend("",  ncol = 2)) + 
# theme(legend.title = element_blank()) + 
geom_text(aes(y = 0.5, label = paste0("n = ", column_total))) + 
scale_y_continuous(breaks=c(0, 0.5, 1), labels = c("0", "0.5", "1"))
save_plot_simple(p_stacked_ct, "CH_stacked_barplot_IDH2_R140_SRSF2_P95", 3.75, 3)


df_ch_mut_plot = ch_mut_df %>% dplyr::filter(ch_mut %in% c("IDH2 SRSF2 MUT", "IDH2 MUT", "SRSF2 MUT"))
df_ch_mut_plot$manual_annotation_v3_l1_final = droplevels(df_ch_mut_plot$manual_annotation_v3_l1_final)
p_stacked_ct = df_ch_mut_plot %>% group_by(donor.id) %>% count(ch_mut) %>% 
  mutate(prop = n / sum(n)) %>%
  mutate(column_total = sum(n)) %>% ungroup() %>% 
ggplot( aes(fill=ch_mut,  x=donor.id, y = prop)) + 
    geom_col(position = position_fill(reverse = TRUE)
    # stat="count"
    ) + scale_fill_manual(values = ch_mut_pal) + #, na.value = makeTransparent("grey"))  +
    ylab("Fraction of Mutant Cells") + xlab('Donor')+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position='bottom')  + 
guides(fill=guide_legend("",  ncol = 2)) + 
# theme(legend.title = element_blank()) + 
geom_text(aes(y = 1.1, label = paste0("n = ", column_total))) + 
scale_y_continuous(breaks=c(0, 0.5, 1), labels = c("0", "0.5", "1"))
save_plot(p_stacked_ct, "CH_stacked_barplot_IDH2_R140_SRSF2_P95_simple", 2, 2)



p_stacked_ct_ch = ggplot(df_ch_mut_plot, aes(fill=manual_annotation_v3_l1_final,  x=donor.id)) + 
    geom_bar(position="fill", stat="count") + scale_fill_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    ylab("Fraction") + xlab('Individual')+ theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
theme(legend.position='bottom') + labs(fill=' ') + 
ggtitle("Cell Type")
# ct_legend = cowplot::get_legend(p_stacked_ct_normal + theme(legend.position = "right")+ guides(fill = guide_legend(ncol = 2, title = "Cell Type")))
save_plot_simple(p_stacked_ct_ch, "stacked_barplot_CH_pre_CRv5_LD_manual_annotation_v3_l1_final", 1.5, 3)



#### CH DGE

DefaultAssay(CH_normal_obj_mono) = "RNA"
meta_de_CH = CH_normal_obj_mono[[]]


library(msigdbr)
library(fgsea)
hallmark_pathway <- gmtPathways("/home/seurat/source_files/h.all.v7.4.symbols.gmt.txt")
head(names(hallmark_pathway))

CD14_mono_CH = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("CH-01", "CH-05")))  %>%  dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono")) %>% tibble::column_to_rownames('cb')
CD14_mono_normals = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("healthy_MSK_02", "healthy_MSK_01"))) %>%   dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono"))  %>% tibble::column_to_rownames('cb')
# HSCs_normals_names = gsub("MSK_", "", rownames(HSCs_normals))
CD14_mono_CH_v_normal = FindMarkers(CH_normal_obj_mono, ident.1 = rownames(CD14_mono_CH), ident.2 = rownames(CD14_mono_normals), logfc.threshold = 0)

CD14_mono_CH_v_normal %>% arrange(avg_logFC)
CD14_mono_CH_v_normal['gene'] = rownames(CD14_mono_CH_v_normal)
write.csv(CD14_mono_CH_v_normal, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH_v_normal_logFC0.csv"), quote = FALSE)
CD14_mono_CH_v_normal = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH_v_normal_logFC0.csv"), row.names = 1)
CD14_mono_CH_v_normal$gene = rownames(CD14_mono_CH_v_normal)
volcano_selected_genes(CD14_mono_CH_v_normal, "CD14_mono_CH_v_normal")


CD14_mono_CH = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("CH-01")))  %>%  dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono")) %>% tibble::column_to_rownames('cb')
CD14_mono_normals = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("healthy_MSK_02", "healthy_MSK_01"))) %>%   dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono"))  %>% tibble::column_to_rownames('cb')
# HSCs_normals_names = gsub("MSK_", "", rownames(HSCs_normals))
CD14_mono_CH_v_normal = FindMarkers(CH_normal_obj_mono, ident.1 = rownames(CD14_mono_CH), ident.2 = rownames(CD14_mono_normals), logfc.threshold = 0)

CD14_mono_CH_v_normal %>% arrange(avg_logFC)
CD14_mono_CH_v_normal['gene'] = rownames(CD14_mono_CH_v_normal)
write.csv(CD14_mono_CH_v_normal, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-01_v_normal_logFC0.csv"), quote = FALSE)
CD14_mono_CH_v_normal = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-01_v_normal_logFC0.csv"), row.names = 1)



CD14_mono_CH = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("CH-05")))  %>%  dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono")) %>% tibble::column_to_rownames('cb')
CD14_mono_normals = meta_de_CH %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( donor.id %in% c("healthy_MSK_02", "healthy_MSK_01"))) %>%   dplyr::filter(manual_annotation_v3_l2_final %in% c("CD14+ Mono"))  %>% tibble::column_to_rownames('cb')
# HSCs_normals_names = gsub("MSK_", "", rownames(HSCs_normals))
CD14_mono_CH_v_normal = FindMarkers(CH_normal_obj_mono, ident.1 = rownames(CD14_mono_CH), ident.2 = rownames(CD14_mono_normals), logfc.threshold = 0)

CD14_mono_CH_v_normal %>% arrange(avg_logFC)
CD14_mono_CH_v_normal['gene'] = rownames(CD14_mono_CH_v_normal)
write.csv(CD14_mono_CH_v_normal, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-05_v_normal_logFC0.csv"), quote = FALSE)
CD14_mono_CH_v_normal = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-05_v_normal_logFC0.csv"), row.names = 1)

outdir = outdir_all
CD14_mono_CH_01_v_normal = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-01_v_normal_logFC0.csv"), row.names = 1)
volcano_selected_genes(CD14_mono_CH_01_v_normal, "CD14_mono_CH_01_v_normal")


CD14_mono_CH_05_v_normal = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_CD14_mono_CH-05_v_normal_logFC0.csv"), row.names = 1)
volcano_selected_genes(CD14_mono_CH_05_v_normal, "CD14_mono_CH_05_v_normal")


### CH DGE 
CD14_mono_CH_v_normal = read.csv(paste0(outdir_all, "CH_vs_healthy_CD14Mono.csv"), row.names = 1)


outdir = outdir_all 

fgsea_res_CH = doFGSEA(CD14_mono_CH_v_normal, "CD14_mono_CH_v_normal", "CH Mono vs Normal Mono")

fgsea_res_CH_plot = fgsea_res_CH %>% dplyr::filter(pathway %in% c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS", "HALLMARK_COMPLEMENT", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_HYPOXIA", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_TGF_BETA_SIGNALING")) %>% 
  mutate(pathway = gsub("_", " ", pathway)) %>% mutate(pathway = gsub("HALLMARK", "", pathway))
p_gsea <- ggplot(fgsea_res_CH_plot, aes(x = NES, y = fct_reorder(pathway, NES))) + 
               geom_point(aes(size = abs(NES), color = padj)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient(limits=c(0, 0.05), low="red") +
        xlab("Normalized Enrichment Score") + 
        ylab("Hallmark Pathway") + # facet_grid(.~type, scales="free") +
        ggtitle("Pathways Enriched in CH Mono")

save_plot_simple(p_gsea, "gsea_dotplot_CH_mono_subset_fromSL", 7, 4)

write.table(fgsea_res_CH, paste0(outdir_all, "CH_GSEA_results.tsv"))
