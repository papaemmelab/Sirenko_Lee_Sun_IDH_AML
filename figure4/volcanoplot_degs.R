#######  volcanoplot 
library(ggrepel)
library(Matrix)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(Seurat)
options(ggrepel.max.overlaps = Inf)

output_dir = "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/integrated/MSK_collaboration/diagnosis/DEGs/"
hallmark_gs <- read.delim("/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/Signatures/h.all.v7.5.1.symbols.gmt", header = F)
hallmark_genesets <- apply(hallmark_gs[,-c(1:2)], MARGIN = 1, function(x){unname(c(x))[x!=""]})
names(hallmark_genesets) <- hallmark_gs[,1]

lineage_list <- readRDS(file = "/gpfs/data/aifantislab/home/sl7424/AML_IDH_NEW/integrated/MSK_collaboration/diagnosis/DEGs/lineage_markers.rds")

# DEG_healthy_vs_AML <- read.csv(paste0(output_dir, "/updated_downsample/DGE_MPP_blasts_v_normal_logFC0_combined_MSK_NYU_downssampled250xpatient_granular_celltype.csv"), row.names = 1)
DEG_NRAS <- read.csv(paste0(output_dir, "/DGE_NRAS_clone_v_parent_clone_MPP_logFC0.csv"), row.names = 1)
# DEG_SRSF2 <- read.csv(paste0(output_dir, "/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0.csv"), row.names = 1)
DEG_SRSF2 <- read.csv(paste0(output_dir, "/DGE_SRSF2_clone_v_parent_clone_MPP_logFC0_srsf2only8_idhonly8.csv"), row.names = 1)



deg_filtered <- DEG_SRSF2

# deg_filtered$Significant <- ifelse(deg_filtered$p_val_adj < 1e-10, "FDR < 1e-10", "Not Sig")
deg_filtered$Significant <- ifelse(deg_filtered$p_val_adj < 1e-5, "FDR < 1e-5", "Not Sig")
deg_filtered$gene <- rownames(deg_filtered)

# red blue pal 
pal_use_man = 
  c("#003366", "#005F73", "#0A9396", # HSC MPP GMP hsc alt  pink alt AD7A99 hsc 2A1C72
    "#BB3E03", "#CA6702", "#F7B2AD", ## ery 
    "#585d86", "#8D91C7", "#866577","#70A9A1", "#77AADD", # myeloid #cd16 champagne D5896F 694F5D #cd14 94D2BD progdc 6E75A8
    "#EE8866", "#D4A6C8", "#99621E", "#805E73" , #B 
    "#555c45" ,  "#828C51" ,  "#A6C36F" ,  "#c99e87" , #NK T pdc A6C36F 335145
    "#005F73")  # eh
names(pal_use_man) = 
  c("HSC",   "MPP" ,  "GMP" ,
    "MEP" , "Ery","Platelet",
    "Progenitor_DC" ,  "cDC", "CD16+ Mono" , "CD14+ Mono" , "Macrophage",  
    "Progenitor_B" ,  "Plasma"  , "B" , "Stroma" ,  
    "CD4+ T" ,"CD8+ T",  "NK",  "pDC", 
    "Immature Mono") 

cellcycle <- c(cc.genes$s.genes, "MKI67")
hypoxia_genes_vol = c("HIF1A", "JUN","STAT3", "STAT5", "IRF1",  "ETS1") #  "FOS", 
UPR_genes_vol = c("XBP1", "ATF3", "DDIT3", "ATF4", "CALR", "ATF6")
nfkb = c("RELA", "RELB", "REL", "RELC", "NFKB2")

myeloid = c(  "SPI1", "CEBPA", "CEBPB","CNRIP1", "FCER1A", "S100A10", "S100A9") # "S100A8"
ery = c( "NFIA","HBD", "HBB", "GATA1", "GATA2")
lymphoid = c("AFF3", "IGLL1", "MZB1", "DNTT", "JCHAIN", "POU2F2", "RAG1", "RAG2")
lineage_HSC = c("CRHBP", "AVP", "PCDH9", "HLF" )
dna_damage = c("MYC", "ATM", "ATR", "RAD50", "RAD51", "RAD52", "E2F1", "CHEK2", "CDK4", "CDK6", "BRCA1", "MDM2", "E2F1")


deg_filtered <- deg_filtered %>% mutate(geneset = ifelse(gene %in% c(hypoxia_genes_vol), "Hypoxia", 
                                                                                     ifelse(gene %in% UPR_genes_vol, "UPR",
                                                                                                   ifelse(gene %in% c(hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION), "OxPHOS", 
                                                                                                   ifelse(gene %in% hallmark_genesets$HALLMARK_TNFA_SIGNALING_VIA_NFKB, "TNFa/NF-kB", 
                                                                                                          ifelse(gene %in% cellcycle, "Cell cycle", NA))))))


geneset_pal = c("#552566","#C29D4E", "darkseagreen4", "#E36414", "goldenrod4", "#1770B8")
names(geneset_pal) <- c("Hypoxia", "UPR", "TNFa/NF-kB", "Cell cycle", "DNA damage",  "OxPHOS")



deg_filtered <- deg_filtered %>% mutate(geneset = ifelse(gene %in% c(hypoxia_genes_vol), "Hypoxia", 
                                                         ifelse(gene %in% lineage_HSC, "HSC", 
                                                                ifelse(gene %in% ery, "Erythroid",
                                                                              ifelse(gene %in% c(hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION), "OxPHOS", 
                                                                                     ifelse(gene %in%  hallmark_genesets$HALLMARK_FATTY_ACID_METABOLISM, "FAM",
                                                                       ifelse(gene %in% c(hallmark_genesets$HALLMARK_TNFA_SIGNALING_VIA_NFKB, "CXCL8"), "TNFa/NF-kB", 
                                                                              ifelse(gene %in% cellcycle, "Cell cycle", 
                                                                                     ifelse(gene %in%  c("CD74", "HOPX","HOXB7", "HOXA5", "HOXA3", "HOXA9"), "Others", NA)))))))))

geneset_pal = c("#552566", "#191970","#C29D4E", "darkseagreen4", "#E36414", "goldenrod4", "#BB3E03", "#1770B8", "black")
names(geneset_pal) <- c("Hypoxia", "HSC" , "UPR", "TNFa/NF-kB", "Cell cycle", "DNA damage",  "FAM", "OxPHOS", "Others")



deg_highlight <- deg_filtered[!is.na(deg_filtered$geneset),] %>% subset(p_val_adj< 1e-2)

NRAS_genes_to_label <- c( "HIF1A", "IRF1", "STAT3", ## hypoxia
                         "SAT1", "AREG", "EIF1", "TNFAIP2", "TLR2", "FOS",  "BTG2",## TNF 
                        "IL18",  "PNRC1", "JUNB", "MAP3K8", "CCNL1", "IFNGR2", "CD44", "CD69",
                        "CEBPD", "ICAM1", "PHLDA1","PLAUR", "DUSP4", ## TNF but high in NRAS
                        "MCM4", "MCM5",  "MCM2", "PCNA", "CDCA7", "GINS2", ## cell cycle
                        "ATF3", "MYC", "HOPX", "CD74", "HOXB7", "HOXA5", "HOXA3", "HOXA9")

oxphos_genes <- hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION[grep(hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION, pattern = "NDUFA|NDUFB|COX")]
oxphos_genes <- intersect(deg_highlight$gene,oxphos_genes)

deg_highlight$gene <- rownames(deg_highlight)
deg_highlight_lable <- subset(deg_highlight, gene %in% c(NRAS_genes_to_label, oxphos_genes, FAM_genes))


ggplot(deg_filtered, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(color = "lightgrey", alpha = 0.8, size = 0.5) +
  theme_classic() + 
  ## gene sets 
  geom_point(data =deg_highlight,  aes(x = avg_log2FC, y = -log10(p_val_adj), color = geneset), , size = 1) +
  theme_classic() + 
  # theme(legend.position = "bottom") +
  geom_text_repel(
    data =  deg_highlight, 
    aes(label = gene, color = geneset),
    size = 2,
    nudge_y   = 0.3,
    nudge_x   = 0.3,
    hjust  = 0, 
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.size = 0.2,
    show.legend=F
  )  + 
  scale_color_manual(values = geneset_pal) 

ggsave(paste0(output_dir, "/", "plot/Volcano_NRAS_parent_MPP_updated.pdf"),  width =6, height =4, dpi = 300, device = "pdf")



### SRSF2

deg_filtered <- DEG_SRSF2

# deg_filtered$Significant <- ifelse(deg_filtered$p_val_adj < 1e-10, "FDR < 1e-10", "Not Sig")
deg_filtered$Significant <- ifelse(deg_filtered$p_val_adj < 0.05, "FDR < 0.05", "Not Sig")
deg_filtered$gene <- rownames(deg_filtered)

APC_genes = rownames(deg_filtered)[grep(rownames(deg_filtered), pattern = "^HLA-")]
deg_filtered <- deg_filtered %>% mutate(geneset = ifelse(gene %in% c(hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION), "OxPHOS",
                                                                              ifelse(gene %in% c(APC_genes, "B2M", "CD74"), "MHC Class I/II",
                                                                                     ifelse(gene %in% c(hallmark_genesets$HALLMARK_TNFA_SIGNALING_VIA_NFKB, "CXCL8"), "TNFa/NF-kB", 
                                                                                            ifelse(gene %in% cellcycle, "Cell cycle",
                                                                                                   ifelse(gene %in% c("CDH9", "FTH1"), "Others",NA))))))

geneset_pal = c("#552566", "#191970","#C29D4E", "darkseagreen4", "#E36414", "#BB3E03", "#1770B8", "black")
names(geneset_pal) <- c("Hypoxia", "HSC" , "MYC", "TNFa/NF-kB", "Cell cycle",  "MHC Class I/II", "OxPHOS", "Others")



deg_highlight <- deg_filtered[!is.na(deg_filtered$geneset),] %>% subset(p_val_adj< 0.05)

SRSF2_genes_to_label <- c("CXCL8", "NR4A1", "NR4A2", "NR4A3" ,"CXCL2", "JUNB", "CEBPB", "AREG", "EIF1",
                          "CD44", "IL18", "MYC", "IL1B",
                          "CDCA7", "TYMS", ## cell cycle
                          "HLA-DPA1", "HLA-DPB1", "CD74", "HLA-DRA", "HLA-DRB5", "HLA-DRB1", "B2M", 
                          "HLA-A", "HLA-B", "HLA-C", 
                          "CDH9", "FTH1")

oxphos_genes <- hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION[grep(hallmark_genesets$HALLMARK_OXIDATIVE_PHOSPHORYLATION, pattern = "NDUFA|NDUFB|COX|LDHA|LDHB")]

deg_highlight$gene <- rownames(deg_highlight)
deg_highlight_label <- subset(deg_highlight, gene %in% c(SRSF2_genes_to_label, oxphos_genes))


ggplot(deg_filtered, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(color = "lightgrey", alpha = 0.8, size = 0.5) +
  theme_classic() + 
  ## gene sets 
  geom_point(data =deg_highlight,  aes(x = avg_log2FC, y = -log10(p_val_adj), color = geneset), size = 1) +
  theme_classic() + 
  # theme(legend.position = "bottom") +
  geom_text_repel(
    data =  deg_highlight_label, 
    aes(label = gene, color = geneset),
    size = 2,
    nudge_y   = 0.3,
    nudge_x   = 0.3,
    hjust  = 0, 
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.size = 0.2,
    show.legend=F
  )  + 
  scale_color_manual(values = geneset_pal) 

# ggsave(paste0(output_dir, "/", "plot/Volcano_SRSF2_parent_MPP_updated.pdf"),  width =6, height =4, dpi = 300, device = "pdf")

ggsave(paste0(output_dir, "/", "plot/Volcano_SRSF2_IDH2_only_gain8.pdf"),  width =6, height =4, dpi = 300, device = "pdf")


