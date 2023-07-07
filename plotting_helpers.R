### plotting helper functions 
library(ggplot2)
library(viridis)
library(ggthemes)
library("RColorBrewer")
library(patchwork)
library(ggpubr)
library(ggridges)
library(EnhancedVolcano)
library(purrr)
library(patchwork)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(scales)




figh = 600
figw = 600
outdir_all = "/home/aifantis_collab/"
fig_outdir = paste0(outdir_all, "manuscript_code/")



### plotting functions 

save_plot = function(p_plot, name_plot, wdth, hght){
  ggsave(plot = p_plot + theme(legend.position='none') , filename = paste0(fig_outdir, name_plot, "_nolegend.png", sep = ""), width = wdth, height = hght, dpi = 500)
  ggsave(plot = p_plot + theme(legend.position='none') , filename = paste0(fig_outdir, name_plot, "_nolegend.pdf", sep = ""), width = wdth, height = hght, dpi = 500)

  ggsave(plot = p_plot + theme(legend.position='bottom') , filename = paste0(fig_outdir, name_plot, "_leg_bottom.png", sep = ""), width = wdth, height = hght, dpi = 500)
  ggsave(plot = p_plot + theme(legend.position='bottom') , filename = paste0(fig_outdir, name_plot, "_leg_bottom.pdf", sep = ""), width = wdth, height = hght, dpi = 500)

  ggsave(plot = p_plot + theme(legend.position='left') , filename = paste0(fig_outdir, name_plot, "_leg_left.png", sep = ""), width = wdth, height = hght, dpi = 500)
  ggsave(plot = p_plot + theme(legend.position='left') , filename = paste0(fig_outdir, name_plot, "_leg_left.pdf", sep = ""), width = wdth, height = hght, dpi = 500)
}

save_plot_simple = function(p_plot, name_plot, wdth, hght){
  ggsave(plot = p_plot , filename = paste0(fig_outdir, name_plot, ".png", sep = ""), width = wdth, height = hght, dpi = 500)
  ggsave(plot = p_plot , filename = paste0(fig_outdir, name_plot, ".pdf", sep = ""), width = wdth, height = hght, dpi = 500)
}


theme_umap = function(){
  theme_minimal() %+replace%  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())     
}
### COLOR PALETTES 

makeTransparent <- function(someColor, alpha=30) scales::alpha(someColor, alpha/100) 

# Cell type palette 

pal_use_man = 
c("#003366", "#005F73", "#0A9396", # HSC MPP GMP 
"#BB3E03", "#CA6702", "#F7B2AD", ## ery 
"#585d86", "#8D91C7", "#866577","#70A9A1", "#77AADD", # mono/dc 
"#EE8866", "#D4A6C8", "#99621E", "#805E73" , #B 
 "#555c45" ,  "#828C51" ,  "#A6C36F" ,  "#c99e87" ) #NK T pdc 
 names(pal_use_man) = 
 c("HSC",   "MPP" ,  "GMP" ,
  "MEP" , "Ery","Platelet",
 "Progenitor_DC" ,  "cDC", "CD16+ Mono" , "CD14+ Mono" , "Macrophage",  
 "Progenitor_B" ,  "Plasma"  , "B" , "Stroma" ,  
 "CD4+ T" ,"CD8+ T",  "NK",  "pDC") 


# donor id color palette
n <- 18
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
patient_pal_v2 <- sample(col_vec, n)
names(patient_pal_v2) = unique(subset(meta_all, (institute_disease %in% c("NYU_AML", "MSK_AML", "MSK_CH")))$donor.id)


### DGE FUNCTIONS 


options(ggrepel.max.overlaps = Inf)
hallmark_pathway <- gmtPathways("/home/seurat/source_files/h.all.v7.4.symbols.gmt.txt")
head(names(hallmark_pathway))



doFGSEA = function(FindMarkers_output, file_name, title){
  # order list, pull out gene name and log2fc, and convert genes to uppercase
  if("avg_log2FC" %in% colnames(FindMarkers_output)) {FindMarkers_output$avg_logFC = FindMarkers_output$avg_log2FC}

  FindMarkers_output <- FindMarkers_output[order(FindMarkers_output$avg_logFC, decreasing = T),]
  FindMarkers_output$Gene.name <- rownames(FindMarkers_output)
  FindMarkers_output <- FindMarkers_output[,c("Gene.name", "avg_logFC")]
  rownames(FindMarkers_output) <- NULL

  FindMarkers_output <- prepare_ranked_list(FindMarkers_output)

  # generate GSEA result table using fgsea() by inputting the pathway list and ranked list
  fgsea_results <- fgsea(pathways = hallmark_pathway,
                    stats = FindMarkers_output,
                    minSize = 15,
                    maxSize = 500 ,
                    nPermSimple = 100000
                    # nperm= 1000
                    )

  p = plot_enrichment(hallmark_pathway, "HALLMARK_HYPOXIA" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir,  "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_HYPOXIA.png"), plot = p)

  p = plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_TNFA_SIGNALING_VIA_NFKB.png"), plot = p)

  p = plot_enrichment(hallmark_pathway, "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_UNFOLDED_PROTEIN_RESPONSE.png"), plot = p)

  p = plot_enrichment(hallmark_pathway, "HALLMARK_MTORC1_SIGNALING" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_MTORC1_SIGNALING.png"), plot = p)

  p = plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_INTERFERON_GAMMA_RESPONSE.png"), plot = p)

  p = plot_enrichment(hallmark_pathway, "HALLMARK_IL6_JAK_STAT3_SIGNALING" , FindMarkers_output)
  ggsave(width=5, height=4,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/enrichment_plot_", file_name, "_DEGS_HALLMARK_IL6_JAK_STAT3_SIGNALING.png"), plot = p)


  p = waterfall_plot(fgsea_results, paste0("Pathways enriched in ", title))
  ggsave(width=6, height=7,filename = paste0(outdir, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/waterfall_", file_name, "_DEGS.png"), plot = p)

  return(fgsea_results)
}


# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_logFC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}
plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}
waterfall_pal = c("#985974", "#679998")
names(waterfall_pal) = c(TRUE, FALSE)
waterfall_plot <- function (fsgea_results__, graph_title) {
  fsgea_results__ %>%
    mutate(short_name = stringr::str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+ scale_fill_manual(values = waterfall_pal) + 
      coord_flip()+ theme_minimal() + 
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7),
            plot.title = element_text(hjust = 1))
}



quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


paletteLength <- 50
myColor <- colorRampPalette(c(muted("blue"), "white", muted("red")))(paletteLength)


plot_heatmap_fsea_results = function(fg_res, celltypes, name){
  fg_prep = fg_res # fg is a list of fgsea results / dataframes 
  # give columns suffix for celltype 
  for(i in cell_types_oi) {
    new_colnames = paste0(colnames(fg_prep[[i]]), "_", i)
    new_colnames = new_colnames[-1]
    colnames(fg_prep[[i]]) = c("pathway", new_colnames)
  }
  fg_all = fg_prep %>% reduce(full_join, by = "pathway") # merge all fg 

  # select only NES columns 
  fg_plot = fg_all %>% select(pathway, starts_with("NES"))
  fg_plot = as.data.frame(fg_plot)
  rownames(fg_plot ) = fg_plot$pathway

  # select specific pathways 
  selected_pathways = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_E2F_TARGETS",  "HALLMARK_G2M_CHECKPOINT" ,  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_P53_PATHWAY", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
  fg_selected = fg_plot %>% dplyr::filter(pathway %in% selected_pathways) 

  # set breaks / scale color 
  to_plot = fg_selected %>% arrange(-NES_MPP) %>% select(-pathway)
  myBreaks <- c(seq(min(fg_selected$NES_MPP), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(fg_selected$NES_MPP)/paletteLength, max(fg_selected$NES_MPP), length.out=floor(paletteLength/2)))
  colnames(to_plot) = gsub("NES_", "",   colnames(to_plot) )   
  colnames(to_plot) = gsub("\\+[[:space:]]", "",   colnames(to_plot) )   # fix name of CD14+ Mono so can use in color annotation below 
  rownames(to_plot) = gsub("HALLMARK_", "", rownames(to_plot))
  rownames(to_plot) = gsub("_", " ", rownames(to_plot))

  # add column annot 
  my_sample_col <- data.frame(celltype = colnames(to_plot))
  row.names(my_sample_col) <- colnames(to_plot)       
  my_colour = list(
      # celltype = c(HSC_MPP = "#5977ff", GMP = "#e17c7c" , CD14Mono = "#D37295", MEP = "#990000")
      #celltype = c(HSC_MPP = "#005F73", GMP = "#0A9396" , CD14Mono = "#70A9A1", MEP_Ery = "#BB3E03")
      celltype = c( MPP = "#005F73", GMP = "#0A9396" , CD14Mono = "#70A9A1", MEP = "#BB3E03") #, Ery = "#CA6702", cDC = "866577", CD16Mono = "#866577", HSC = "#003366",)

  )

  ## get pvalues for annotation 
  fg_pval = fg_all %>% select(pathway, starts_with("padj"))
  fg_pval = as.data.frame(fg_pval)
  rownames(fg_pval ) = fg_pval$pathway
  fg_pval_selected = fg_pval %>% dplyr::filter(pathway %in% selected_pathways) %>% select(-pathway)
  colnames(fg_pval_selected) = gsub("padj_", "",   colnames(fg_pval_selected) )   
  colnames(fg_pval_selected) = gsub("\\+[[:space:]]", "",   colnames(fg_pval_selected) )   # fix name of CD14+ Mono so can use in color annotation below 
  rownames(fg_pval_selected) = gsub("HALLMARK_", "", rownames(fg_pval_selected))
  rownames(fg_pval_selected) = gsub("_", " ", rownames(fg_pval_selected))


  fg_pval_sorted = fg_pval_selected[match(rownames(to_plot), rownames(fg_pval_selected)), ]
  fg_pval_sorted[fg_pval_sorted<0.05] = "*"
  fg_pval_sorted[fg_pval_sorted>0.05] = ""

  #plot 
  png(width = figw/1.5, height = figh/2, pointsize = 1,  filename = paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/", "heatmap_fgsea_", name, "_selected.png", sep = ""))
  p = pheatmap(to_plot, cluster_cols = F, cluster_rows = F, breaks = myBreaks, color = myColor, annotation_col = my_sample_col, annotation_colors = my_colour,  display_numbers = fg_pval_sorted, number_color = "white")
  dev.off()

  pdf(width = 5, height = 4, pointsize = 1,  file = paste0(outdir_all, "de_outdir_manual_annot_CRv5_seurat_int_man_annot/de_figures/", "heatmap_fgsea_", name, "_selected.pdf", sep = ""))
  p = pheatmap(to_plot, cluster_cols = F, cluster_rows = F, breaks = myBreaks, color = myColor, annotation_col = my_sample_col, annotation_colors = my_colour,  display_numbers = fg_pval_sorted, number_color = "white")
  dev.off()

}



