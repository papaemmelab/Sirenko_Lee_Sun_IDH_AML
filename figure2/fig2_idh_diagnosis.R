
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
source("../figure1/fig1_umaps.R")


##### 
meta_dx_norm = meta_all0 %>% dplyr::filter(manual_annotation_v3_l1_final %in% c("CD14+ Mono","CD16+ Mono","HSC","MPP",
"GMP","MEP","cDC"  , "Progenitor_DC" ,"Ery", "Macrophage"     ))  %>% dplyr::filter(stage_simple %in% c("healthy", "diagnosis"))


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
meta_dx_norm = meta_dx_norm %>% mutate(donor.id.original = donor.id) 
for (i in 1:length(names_dict)){meta_dx_norm$donor.id <- replace(meta_dx_norm$donor.id, meta_dx_norm$donor.id == names(names_dict[i]), as.character(names_dict[i]))}


# UMAP OF NORMALs 

p_all_umap = ggplot(meta_dx_norm,  aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point(data = subset(meta_dx_norm, (institute_disease %in% c("MSK_healthy", "NYU_healthy")))),  aes(color = as.factor(manual_annotation_v3_l1_final) , size=1.5) + 
    scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    theme(legend.position='bottom') + 
    ggtitle("Normal Ref UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = 2)))
save_plot(p_all_umap, "umap_normal_ref", 8, 8)


# COLOR DIAGNOSIS UMAP BY CELLTYPE 

p_all_umap = ggplot(meta_dx_norm,  aes(x=rpca_UMAP_1, y=rpca_UMAP_2) ) + theme_minimal() + #xlim(-10, 12) + ylim(-7.5, 5) + 
    geom_point(data = subset(meta_dx_norm, (institute_disease %in% c("MSK_healthy", "NYU_healthy"))), colour = makeTransparent("grey") , size=1.5) + 
    geom_point(data = subset(meta_dx_norm, (institute_disease %in% c("MSK_AML", "NYU_AML"))),  aes(color = as.factor(manual_annotation_v3_l1_final)) , #pch=21,
    size=umap_pt_size)+ scale_color_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    theme_umap() + 
    theme(legend.position='bottom') + 
    ggtitle("Diagnosis UMAP")+ ylab("UMAP 2") + xlab('UMAP 1') + labs(color=' ') + guides(color = guide_legend(override.aes = list(size = 2)))
save_plot(p_all_umap, "umap_AML_dx", 8, 8)



### MAKE STACKED BARPLOT FOR CELLTYPE REPRESENTATION 

meta_dx_norm = meta_dx_norm %>% mutate(manual_annotation_v3_l1_final = factor(manual_annotation_v3_l1_final, levels = c("HSC", "MPP",  "GMP", "CD14+ Mono", "CD16+ Mono", "Macrophage", "Progenitor_DC", "cDC", "MEP", "Ery")))
p_stacked_ct_aml = ggplot(meta_dx_norm, aes(fill=forcats::fct_rev(manual_annotation_v3_l1_final),  x=donor.id)) + 
    geom_bar(position="fill", stat="count") + scale_fill_manual(values = pal_use_man, na.value = makeTransparent("grey"))  +
    ylab("Fraction") + xlab('Individual')+ theme_bw() + facet_grid(. ~ group, scales = "free", space='free')  + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
theme(legend.position='bottom') + labs(fill=' ') + 
ggtitle("Cell Type")
save_plot_simple(p_stacked_ct_aml, "stacked_barplot_blast_subtypes_AML_pre_CRv5_LD_manual_annotation_v3_l1_final", 8, 4)


### BOXPLOT FRACTION IMMATURE CELLS 

df_mye_ratio_total = meta_dx_norm %>% mutate(myeloid = case_when(manual_annotation_v3_l1_final %in% c("CD14+ Mono","CD16+ Mono", "cDC" ,"Macrophage"  ) ~ "mature",
                                                            manual_annotation_v3_l1_final %in% c("HSC","MPP","GMP","MEP", "Progenitor_DC", "Ery" ) ~ "immature")) %>% 
                                                             dplyr::count(donor.id.original, myeloid) %>% group_by (donor.id.original) %>% summarise(mye_ratio_total = n[myeloid=="immature"]/sum(n))
p_box_mye = df_mye_ratio_total %>% mutate(disease = case_when(donor.id.original %in% c("BM_healthy_NYU_1" ,"BM_healthy_NYU_2" ,"BM_healthy_NYU_3",
 "BM_healthy_NYU_4", "BM_healthy_NYU_5" ,"BM_healthy_NYU_6" ,"BM_healthy_NYU_7","BM_healthy_NYU_8", "healthy_MSK_01"  , "healthy_MSK_02" ,  "healthy_MSK_03" ,"healthy_MSK_04" ) ~ "healthy",
                                            donor.id.original %in% c("AML0024","AML0048","AML0875","AML1371", "AML1529","AML1896","AML2123","AML4340","AML4897","IDH1i-01","IDH1i-02", "IDH1i-03"  ,      
 "IDH2i-01","IDH2i-02","IDH2i-03"  ) ~ "AML"))  %>% mutate(disease = factor(disease, levels = c("healthy", "AML"))) %>% 
                                            ggplot(aes(x = disease, y = mye_ratio_total)) + geom_boxplot(outlier.shape = NA, fill = "lightblue") + geom_jitter(color="black", size=0.4, alpha=0.9) +
                                            stat_compare_means(size = 2) + ylim(0,1.2) + 
                                            theme_classic() + scale_x_discrete(breaks=c( "AML", "healthy"),labels=c("AML", "Healthy")) + ylab("Fraction Immature / Total Myeloid") + xlab("") + 
                                            theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 save_plot_simple(p_box_mye, "boxplot_myeloid_ratio_dx_total_stat", 2, 3)                                          

### SCORE CELLS FOR EXPRESSION OF STEMNESS AND LINEAGE MODULES 


##### add in van galen module scores 
vanGalen_markers_df = read.csv("/home/seurat/source_files/vanGalen2019_genesets.csv")
vanGalen_modules = colnames(vanGalen_markers_df)

calc_module_score_seurat_from_msigdbr = function(seurat_object_plot, name_seurat_object, name_module){
    module = as.list(msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_name %in% name_module) %>% dplyr::select(human_gene_symbol))

    seurat_object_plot = calc_module_score_seurat(seurat_object_plot, name_seurat_object, module, name_module )
}
calc_module_score_seurat = function(seurat_object_plot, name_seurat_object, module, name_module){
    seurat_object_plot <- AddModuleScore(
    object = seurat_object_plot,
    features = module,
    ctrl = 50,
    nbin = 10,
    name = paste0(name_module, "_seurat.score"), 
    seed = 123, 
    search = TRUE,
    assay = "RNA"
    )
    return(seurat_object_plot)
}


colnames(vanGalen_markers_df) = paste0("vanGalen19_", colnames(vanGalen_markers_df))

for (i in seq_along(colnames(vanGalen_markers_df))){
    module = na.omit(data.frame("module"=unname(unlist(data.frame(t(unique(vanGalen_markers_df[colnames(vanGalen_markers_df)[i]]))))), stringsAsFactors = FALSE))
    module = module[module != "?"]
    module = module[module != ""]
    print(colnames(vanGalen_markers_df)[i])
    print(module)
    # module = as.list(module)
    mod1 = list()
    mod1[[colnames(vanGalen_markers_df)[i]]] = module
    # sampls_seurat_query = modify2( sampls_seurat_query, "sampls_seurat_query", calc_module_score_seurat,  module = mod1, name_module = colnames(vanGalen_markers_df)[i])
    int_obj_all = calc_module_score_seurat(int_obj_all, module = mod1, name_module = colnames(vanGalen_markers_df)[i])
}


Velten_priming = read.csv("/home/seurat/patient_centric_subset/source_files/Velten_supp_table4_lineage_priming.csv")
Velten_priming = Velten_priming %>% mutate_all(str_remove_all,  "\\([^()]+?\\)")
Velten_priming = Velten_priming %>% mutate_all(str_trim, "right")
Velten_priming = Velten_priming %>% mutate_all(str_replace_all,  "\\?.*\\)", "")
Velten_priming = Velten_priming %>% mutate_all(str_remove_all,  "U+003F")


for (i in seq_along(colnames(Velten_priming))){
    module = na.omit(data.frame("module"=unname(unlist(data.frame(t(unique(Velten_priming[colnames(Velten_priming)[i]]))))), stringsAsFactors = FALSE))
    module = module[module != "?"]
    module = module[module != ""]
    print(colnames(Velten_priming)[i])
    print(module)
    # module = as.list(module)
    mod1 = list()
    mod1[[colnames(Velten_priming)[i]]] = module
    # sampls_seurat_v2 = modify2( sampls_seurat_v2, names(sampls_seurat_v2), calc_module_score_seurat,  module = mod1, name_module = colnames(Velten_priming)[i])
    # sampls_seurat_query = modify2( sampls_seurat_query, names(sampls_seurat_query), calc_module_score_customsets,  module = mod1, name_module = colnames(Velten_priming)[i])
    int_obj_all = calc_module_score_seurat(int_obj_all, "int_obj_all", module = mod1, name_module = colnames(Velten_priming)[i])

}

subset_dx_blasts = subset(subset(int_obj_all, stage_simple == "diagnosis"), manual_annotation_v3_l1_final %in% c("Immature Mono", "CD14+ Mono",    "MPP","cDC","GMP", "Progenitor_DC" ,   "Ery", "MEP","HSC","CD16+ Mono","Macrophage"     ) ) 
#subset_dx_blasts = preprocess_seurat(subset_dx_blasts)



df_score_plot = subset_dx_blasts[[]]
df_score_plot = df_score_plot %>% mutate(donor.id.original = donor.id) 
for (i in 1:length(names_dict)){df_score_plot$donor.id <- replace(df_score_plot$donor.id, df_score_plot$donor.id == names(names_dict[i]), as.character(names_dict[i]))}


df_score_plot$TAL1.HSF1_ery_score_scaled <- scales::rescale(df_score_plot$TAL1.HSF1_seurat.score1, to=c(0,1))
df_score_plot$vanGalen19_Normal_HSC_Prog_score_scaled <- scales::rescale(df_score_plot$vanGalen19_Normal_HSC_Prog_seurat.score1, to=c(-0.8,0.5))
df_score_plot$vanGalen19_Normal_GMP_score_scaled <- scales::rescale(df_score_plot$vanGalen19_Normal_GMP_seurat.score1, to=c(0,1))
df_score_plot$vanGalen19_AML_Myeloid.like_score_scaled <- scales::rescale(df_score_plot$vanGalen19_AML_Myeloid.like_seurat.score1, to=c(0,1))

df_score_plot = df_score_plot %>% mutate(mye_ery_split = case_when(TAL1.HSF1_ery_score_scaled >= vanGalen19_Normal_GMP_score_scaled ~ -1* TAL1.HSF1_ery_score_scaled,
                                                                   TAL1.HSF1_ery_score_scaled < vanGalen19_Normal_GMP_score_scaled ~ vanGalen19_Normal_GMP_score_scaled ))


df_score_plot$cell_type = df_score_plot$manual_annotation_v3_l1_final
df_score_plot = df_score_plot %>% dplyr::filter(!(cell_type %in% c("CD14+ Mono", "CD16+ Mono", "cDC", "Progenitor_DC", "Macrophage", "Immature Mono")))

df_score_plot_fork = df_score_plot #%>% dplyr::filter(!(cell_type %in% c("CD14+ Mono", "CD16+ Mono", "cDC", "Progenitor_DC", "Macrophage")))
df_score_plot_fork = df_score_plot_fork %>% mutate(gmp_ery_dif = vanGalen19_Normal_GMP_score_scaled - TAL1.HSF1_ery_score_scaled )

# set pt order to plot 
order_pts = c("IDH1i-01", "IDH1i-03", "AML1529" , "AML4897"  , "AML1371" , "IDH2i-01", "IDH2i-02","IDH2i-03",  "AML0875" , "IDH1i-02", "AML0048", "AML2123", "AML4340",  "AML0024", "AML1896" )
for (i in 1:length(names_dict)){order_pts <- replace(order_pts, order_pts == names(names_dict[i]), as.character(names_dict[i]))}


df_score_plot_fork$donor.id = factor(df_score_plot_fork$donor.id, levels = order_pts)
df_score_plot_fork0 = df_score_plot_fork %>% mutate(jitter = (rnorm(n(), mean=gmp_ery_dif)) /100 ) #-0.15 + (sample(0:300, n(), replace = TRUE, rnorm(n, mean)))/1000)
df_score_plot_fork0  = df_score_plot_fork0 %>% mutate(mye_ery_score_processed = case_when(vanGalen19_Normal_HSC_Prog_score_scaled <= 0 ~ gmp_ery_dif,
                                                                                  vanGalen19_Normal_HSC_Prog_score_scaled > 0 ~ 2*jitter))
df_score_plot_fork0  = df_score_plot_fork0 %>% mutate(stem_score_processed = case_when(vanGalen19_Normal_HSC_Prog_score_scaled > 0 ~ vanGalen19_Normal_HSC_Prog_score_scaled,
                                                                                        gmp_ery_dif <= 0 ~ (gmp_ery_dif) + 2*jitter,
                                                                                        gmp_ery_dif > 0 ~ -(gmp_ery_dif) + 2*jitter))

df_score_plot_fork0$cell_type = factor(df_score_plot_fork0$cell_type, levels = c( "GMP", "Ery", "MEP", "MPP", "HSC"))


p_fork = ggplot(df_score_plot_fork0, aes(x = mye_ery_score_processed, y = stem_score_processed, color = cell_type)) + geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0.5), color = "darkgrey") + geom_segment(aes(x = 0, xend = 1, y = 0, yend = -1), color = "darkgrey") + geom_segment(aes(x = 0, xend = -1, y = 0, yend = -1), color = "darkgrey") + 
geom_point(size = 0.01) + scale_color_manual(values = pal_use_man) + 
theme_void() + ylab("") + xlab("") + theme(legend.position='none') 
save_plot_simple(p_fork, "fork_stemness_v1_all", 1.5,1.5)



NPM1_group = c("IDH1i-01", "IDH1i-03", "AML1529" , "AML4897" ) # , "AML1371" )
SRSF2_group = c("IDH2i-01", "IDH2i-02","IDH2i-03",  "AML0875"  )
IDH1_group = c("IDH1i-01", "IDH1i-02", "IDH1i-03",  "AML2123", "AML4340", "AML1529") # "AML0048",
IDH2_group = c("IDH2i-01", "IDH2i-02","IDH2i-03",  "AML4897", "AML1371", "AML0875", "AML0024", "AML1896")
for (i in 1:length(names_dict)){NPM1_group <- replace(NPM1_group, NPM1_group == names(names_dict[i]), as.character(names_dict[i]))}
for (i in 1:length(names_dict)){SRSF2_group <- replace(SRSF2_group, SRSF2_group == names(names_dict[i]), as.character(names_dict[i]))}
for (i in 1:length(names_dict)){IDH1_group <- replace(IDH1_group, IDH1_group == names(names_dict[i]), as.character(names_dict[i]))}
for (i in 1:length(names_dict)){IDH2_group <- replace(IDH2_group, IDH2_group == names(names_dict[i]), as.character(names_dict[i]))}



df_score_plot_fork0 = df_score_plot_fork0 %>% mutate(genotype_group = case_when(donor.id %in%  NPM1_group ~ "NPM1", 
                                                                    donor.id %in%  SRSF2_group ~ "SRSF2", 
                                                                    TRUE ~ "Other"))


df_score_plot_fork0 = df_score_plot_fork0 %>% mutate(IDH_group = case_when(donor.id %in% IDH1_group ~ "IDH1", 
                                                                    donor.id %in% IDH2_group ~ "IDH2", 
                                                                    TRUE ~ "Other"))



p_fork = ggplot(df_score_plot_fork0, aes(x = mye_ery_score_processed, y = stem_score_processed, color = cell_type)) + geom_segment(aes(x = 0, xend = 0, y = 0, yend = 0.5), color = "darkgrey") + geom_segment(aes(x = 0, xend = 1, y = 0, yend = -1), color = "darkgrey") + geom_segment(aes(x = 0, xend = -1, y = 0, yend = -1), color = "darkgrey") + 
geom_point(size = 0.01) + scale_color_manual(values = pal_use_man) + 
theme_void() + ylab("") + xlab("") + theme(legend.position='none') + facet_wrap(.~genotype_group)
save_plot_simple(p_fork, "fork_stemness_v1_all_split_genotype_group_axes", 4.5,1.5)



p_fork = ggplot(df_score_plot_fork0, aes(x = mye_ery_score_processed, y = stem_score_processed, color = cell_type)) + geom_point(size = 0.01) + scale_color_manual(values = pal_use_man) + 
theme_void() + ylab("") + xlab("") + theme(legend.position='none') + facet_wrap(.~IDH_group)
save_plot_simple(p_fork, "fork_stemness_v1_all_split_IDH", 3,1.5)




## add density ridge plots 
genotype_group_cols = c("#AA7A38",  "#0c4f4e", "lightgray")
names(genotype_group_cols) = c("NPM1", "SRSF2", "Other")
df_ridgeplot = df_score_plot_fork0
df_ridgeplot$cell_type = factor(df_ridgeplot$cell_type, levels = c( "Ery", "MEP", "GMP", "MPP", "HSC"))

p_ridge = df_ridgeplot %>% dplyr::filter(vanGalen19_Normal_HSC_Prog_score_scaled <= 0, genotype_group %in% c("NPM1", "SRSF2")) %>% ggplot( aes(x = mye_ery_score_processed, y = genotype_group, fill = genotype_group)) + geom_density_ridges(color = "white") + 
  theme_classic() + scale_fill_manual(values = genotype_group_cols) + xlab("Myeloid - Ery Score") + ylab("Genotype")+ theme(legend.position='none') + xlim(-1,1)
save_plot_simple(p_ridge, "ridgeplot_mye_ery_score_NPM1_SRSF2", 5,1.5)


p_ridge = df_ridgeplot %>% dplyr::filter(vanGalen19_Normal_HSC_Prog_score_scaled <= 0, IDH_group %in% c("IDH1", "IDH2")) %>% ggplot( aes(x = mye_ery_score_processed, y = IDH_group, fill = IDH_group)) + geom_density_ridges(color = "white") + 
  theme_classic() + scale_fill_manual(values = c( "#b64c63", "#5c315f"))  + xlab("Myeloid - Ery Score") + ylab("Genotype")+ theme(legend.position='none') + xlim(-1,1)
save_plot_simple(p_ridge, "ridgeplot_mye_ery_score_IDH1_IDH2", 5,1.5)


p_ridge = df_ridgeplot  %>% ggplot( aes(x = mye_ery_score_processed, y = cell_type, fill = cell_type)) + geom_density_ridges(color = "white") + 
  theme_classic() + scale_fill_manual(values = pal_use_man) + xlab("Myeloid - Ery Score") + ylab("Cell Type")+ theme(legend.position='none') + xlim(-1,1)
save_plot_simple(p_ridge, "ridgeplot_celltype_mye_ery_nooutline", 5,2)

df_ridgeplot$cell_type = factor(df_ridgeplot$cell_type, levels = c("HSC", "MPP", "GMP", "MEP", "Ery"))
p_ridge = df_ridgeplot  %>% ggplot( aes(x = vanGalen19_Normal_HSC_Prog_score_scaled, y = cell_type, fill = cell_type)) + geom_density_ridges(color = "white") + 
  theme_classic() + scale_fill_manual(values = pal_use_man) + xlab("Stemness") + ylab("Cell Type") + coord_flip()+ theme(legend.position='none') 
save_plot_simple(p_ridge, "ridgeplot_celltype_stemness_nooutline", 1.5, 5)



### BOXPLOT: Lineage Scores in NPM1 vs SRSF2 
p_box_score = df_score_plot_fork0 %>% filter(genotype_group %in% c("NPM1", "SRSF2"), cell_type == "MPP") %>% ggplot( aes(x = genotype_group, y = TAL1.HSF1_ery_score_scaled, fill = genotype_group)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + ylab("MEP/Ery Score") + xlab("Genotype") + scale_fill_manual(values = c( "#AA7A38", "#0c4f4e")) + theme(legend.position='none') + ggtitle("MPP") + scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) + 
  stat_compare_means(label.sep = "\n",  label.y = 0.8, size = 2)
save_plot_simple(p_box_score, "boxplot_ery_score_MPP_scaled", 1.5, 2)

p_box_score = df_score_plot_fork0 %>% filter(genotype_group %in% c("NPM1", "SRSF2"), cell_type == "MPP") %>% ggplot( aes(x = genotype_group, y = vanGalen19_Normal_GMP_score_scaled, fill = genotype_group)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + ylab("GMP Score") + xlab("Genotype") + scale_fill_manual(values = c( "#AA7A38", "#0c4f4e")) + theme(legend.position='none') + ggtitle("MPP") + scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1))+ 
    stat_compare_means(label.sep = "\n",  label.y = 0.8, size = 2)
save_plot_simple(p_box_score, "boxplot_gmp_score_MPP_scaled", 1.5, 2)

test_npm1 = df_score_plot_fork0 %>% filter(genotype_group %in% c("NPM1"), cell_type == "MPP")
test_srsf2 = df_score_plot_fork0 %>% filter(genotype_group %in% c("SRSF2"), cell_type == "MPP")
t.test(test_npm1$TAL1.HSF1_seurat.score1, test_srsf2$TAL1.HSF1_seurat.score1)
t.test(test_npm1$vanGalen19_Normal_GMP_seurat.score1, test_srsf2$vanGalen19_Normal_GMP_seurat.score1)


### BOXPLOT: Lineage Scores in IDH1 vs IDH2 
p_box_score = df_score_plot_fork0 %>% filter(IDH_group %in% c("IDH1", "IDH2"), cell_type == "MPP") %>% ggplot( aes(x = IDH_group, y = TAL1.HSF1_ery_score_scaled, fill = IDH_group)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + ylab("MEP/Ery Score") + xlab("Genotype") + scale_fill_manual(values = c( "#b64c63", "#5c315f")) + theme(legend.position='none') + ggtitle("MPP") + scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) + 
  stat_compare_means(label.sep = "\n",  label.y = 0.8, size = 2)
save_plot_simple(p_box_score, "boxplot_ery_score_MPP_scaled_IDH1_IDH2", 1.5, 2)

p_box_score = df_score_plot_fork0 %>% filter(IDH_group %in% c("IDH1", "IDH2"), cell_type == "MPP") %>% ggplot( aes(x = IDH_group, y = vanGalen19_Normal_GMP_score_scaled, fill = IDH_group)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + ylab("GMP Score") + xlab("Genotype") + scale_fill_manual(values = c( "#b64c63", "#5c315f")) + theme(legend.position='none') + ggtitle("MPP") + scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1))+ 
    stat_compare_means(label.sep = "\n",  label.y = 0.8, size = 2)
save_plot_simple(p_box_score, "boxplot_gmp_score_MPP_scaled_IDH1_IDH2", 1.5, 2)


p_box_score = df_score_plot_fork0 %>% filter(IDH_group %in% c("IDH1", "IDH2"), cell_type == "MPP") %>% ggplot( aes(x = IDH_group, y = vanGalen19_Normal_HSC_Prog_seurat.score1, fill = IDH_group)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + ylab("Stmness Score") + xlab("Genotype") + scale_fill_manual(values = c("#b64c63", "#5c315f")) + theme(legend.position='none') + ggtitle("MPP") + scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1))+ 
    stat_compare_means(label.sep = "\n",  label.y = 0.8, size = 2)
save_plot_simple(p_box_score, "boxplot_HSC_score_MPP_scaled_IDH1_IDH2", 1.5, 2)


test_idh1 = df_score_plot_fork0 %>% filter(IDH_group %in% c("IDH1"), cell_type == "MPP")
test_idh2 = df_score_plot_fork0 %>% filter(IDH_group %in% c("IDH2"), cell_type == "MPP")
t.test(test_idh1$TAL1.HSF1_seurat.score1, test_idh2$TAL1.HSF1_seurat.score1)
t.test(test_idh1$vanGalen19_Normal_GMP_seurat.score1, test_idh2$vanGalen19_Normal_GMP_seurat.score1)


### DIFFERENTIAL GENE EXPRESSION

#### AML vs Normal Comparisons
## Downsample AML vs Normal by patient 
meta_all0_dx_norm_MPP = meta_all0 %>% dplyr::filter(stage_simple %in% c("diagnosis", "healthy"), manual_annotation_v3_l1_final %in% c("MPP"))
dx_MPP_obj = subset(int_obj_all, cells = rownames(meta_all0_dx_norm_MPP)) 

Idents(dx_MPP_obj) = "donor.id"
dx_norm_MPP_small <- subset(dx_MPP_obj, downsample = 500)

Idents(dx_norm_MPP_small) = "stage_simple"
DefaultAssay(dx_norm_MPP_small) = "RNA"
MPP_AML_vs_norm_downsampled = FindMarkers(dx_norm_MPP_small, ident.1 = "diagnosis", ident.2 = "healthy", logfc.threshold = 0.1, assay = "RNA")
MPP_AML_vs_norm_downsampled$gene = rownames(MPP_AML_vs_norm_downsampled)   


### AML vs normal comparison across all cell types 

de_subset_obj = subset(int_obj_all, (stage_simple %in% c("healthy", "diagnosis")))
Idents(de_subset_obj) = "stage_simple"
DefaultAssay(de_subset_obj) = "RNA"

de_subset_obj_meta = de_subset_obj[[]]
de_subset_obj_meta = de_subset_obj_meta %>% mutate(group_DE = case_when(institute_disease %in% c("MSK_AML", "NYU_AML") ~ "AML", 
                                                                                  institute_disease %in% c("NYU_healthy" ,"MSK_healthy") ~ "Healthy",
                                                                                  TRUE ~ "other")) %>% mutate(cell_type_de = case_when(manual_annotation_v3_l1_final %in% c("HSC", "MPP") ~ "HSC_MPP",
                                                                                                                                        manual_annotation_v3_l1_final %in% c("MEP", "Ery") ~ "MEP_Ery",
                                                                                                                                        TRUE ~ as.character(manual_annotation_v3_l1_final)))

de_subset_obj_meta = de_subset_obj_meta %>% mutate(group_DE = case_when(institute_disease %in% c("MSK_AML", "NYU_AML") ~ "AML", 
                                                                                  institute_disease %in% c("NYU_healthy" ,"MSK_healthy") ~ "Healthy",
                                                                                  TRUE ~ "other")) %>% mutate(cell_type_de_granular = case_when(manual_annotation_v3_l1_final %in% c("Immature Mono") ~ "MPP",
                                                                                                                                        TRUE ~ as.character(manual_annotation_v3_l1_final)))



cell_types_oi = c( "MPP", "GMP", "CD14+ Mono", "MEP") # , "Ery", "HSC","CD16+ Mono", "cDC", )#c("MEP_Ery" ) #"HSC_MPP", "GMP", "CD14+ Mono")

for(i in cell_types_oi) {
  blasts_subset = de_subset_obj_meta %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( institute_disease %in% c("MSK_AML", "NYU_AML")))  %>%  dplyr::filter(cell_type_de_granular %in% c(i)) %>% dplyr::filter(stage_simple == "diagnosis")   %>% tibble::column_to_rownames('cb')
  normals_subset = de_subset_obj_meta %>% tibble::rownames_to_column('cb') %>% dplyr::filter(( institute_disease %in% c("NYU_healthy" ,"MSK_healthy"))) %>%   dplyr::filter(cell_type_de_granular %in% c(i))  %>% tibble::column_to_rownames('cb')
  set.seed(111)
  blasts_downsampl = blasts_subset %>% tibble::rownames_to_column('cb') %>% group_by(donor.id) %>% dplyr::slice_sample(n = 500) %>% tibble::column_to_rownames('cb')
  normals_downsampl = normals_subset %>% tibble::rownames_to_column('cb') %>% group_by(donor.id) %>% dplyr::slice_sample(n = 500) %>% tibble::column_to_rownames('cb')
  
  test_obj = subset(de_subset_obj, cells = c(rownames(blasts_downsampl), rownames(normals_downsampl)))
  Idents(test_obj) = "stage_simple"
  DefaultAssay(test_obj) = "RNA"

  markers[[i]] = FindMarkers(test_obj, ident.1 = "diagnosis", ident.2 = "healthy", logfc.threshold = 0, assay = "RNA")
  markers[[i]]['gene'] = rownames(markers[[i]])
  write.csv(markers[[i]], paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_", i, "_blasts_v_normal_logFC0_combined_MSK_NYU_downssampled500xpatient_granular_celltype.csv"), quote = FALSE)

}
saveRDS(markers, paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_marker_list_blasts_v_normal_logFC0_combined_MSK_NYU_downssampled500xpatient_granular_celltype.rds"))
markers = readRDS(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_marker_list_blasts_v_normal_logFC0_combined_MSK_NYU_downssampled500xpatient_granular_celltype.rds"))



### HEATMAP GSEA RESULTS 

markers_read = list()
fg = list()
fg_go = list()
for(i in cell_types_oi) {
  # i = "MPP" 
  markers_read[[i]] = read.csv(paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/DGE_", i, "_blasts_v_normal_logFC0_combined_MSK_NYU_downssampled250xpatient_granular_celltype.csv"), row.names = 1)
  markers_read[[i]]$avg_logFC = markers_read[[i]]$avg_2logFC

  markers_read[[i]] = markers_read[[i]] %>% dplyr::filter(avg_log2FC != 0 ) 
  markers_read[[i]] = markers_read[[i]] %>% dplyr::filter(p_val_adj != -Inf ) 
  markers_read[[i]] = markers_read[[i]] %>% dplyr::filter(avg_log2FC != Inf ) 
  # markers_read[[i]] = markers_read[[i]][Reduce(`&`, lapply(markers_read[[i]], is.finite)),]
  markers_read[[i]] = na.omit(markers_read[[i]])
  fg[[i]] = doFGSEA(markers_read[[i]] , paste0(i, "_AML_vs_normal_MSK_NYU_downsampled250pt"), paste0(i, "_AML_vs_normal_MSK_NYU_downsampled250pt_granular_celltype"))
    fg_go[[i]] = doFGSEA_go(markers_read[[i]] , paste0(i, "_AML_vs_normal_MSK_NYU_downsampled250pt"), paste0(i, "_AML_vs_normal_MSK_NYU_downsampled250pt_granular_celltype"))

}

name = "MPP"
fg_selected = fg[[name]] %>% dplyr::filter(pathway %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_E2F_TARGETS",  "HALLMARK_G2M_CHECKPOINT" ,  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")) 
fg_selected = as.data.frame(fg_selected)
fg_selected = fg_selected %>% select(pathway, NES)
rownames(fg_selected ) = fg_selected$pathway
to_plot = fg_selected %>% select(NES) %>% mutate(NES = -NES)%>% arrange(-NES) 
myBreaks <- c(seq(min(fg_selected$NES), 0, length.out=ceiling(paletteLength/2) + 1), 
            seq(max(fg_selected$NES)/paletteLength, max(fg_selected$NES), length.out=floor(paletteLength/2)))
png(width = figw, height = figh, pointsize = 1,  filename = paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/", "heatmap_fgsea_", name, "_selected_downsampled250pt_granular_celltype.png", sep = ""))
p = pheatmap(to_plot, cluster_cols = F, cluster_rows = F, breaks = myBreaks, color = myColor)
dev.off()


plot_heatmap_fsea_results(fg, cell_types_oi, "combined_downsampl250_granular_celltype")





