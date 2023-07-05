
#### Co-mutation analysis in Tazi et al and Papaemmanuil et al. 

library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(tidyr)


##### Plot the odds ratios 
gene_gene_interactions_OR <- function(genomicData,columns,genomic_groups,column_order=FALSE){
  # tmp = merged_cohorts
  # genomicData = tmp
  # columns=c("IDH1","IDH2_p.R140","IDH2_p.R172", "DNMT3A","SRSF2", "NPM1", "NRAS_p.G12_13", "NRAS_p.Q61_62", "TET2", "RUNX1", "ASXL1", "STAG2", "complex", "add_8",  "KRAS", "ITD", "TP53", "complex", "del_7", "PTPN11")
  # genomic_groups="Genetics"
  # column_order=T
  
  
  genomicData <- genomicData[,columns]
  genomicData <- (sapply(unique(sub("(_ITD)|(_TKD)|(_other)|(_mono)|(_bi)|(_p172)|(_p140)","",colnames(genomicData) )), function(x) rowSums(genomicData[grep(paste(x,"($|_.+)",sep=""), colnames(genomicData))])) > 0)+0
  genomicData <- genomicData[,colSums(genomicData, na.rm=TRUE)>=8]
  
  # update columns if removed from genomicData 
  columns = columns[ columns %in% colnames(genomicData)]
  genomicGroups <- factor(grepl("^[a-z]", colnames(genomicData)) + grepl("t_*[0-9M]", colnames(genomicData) ), labels=genomic_groups)
  dim(genomicData)
  #dev.new(width=20,height=20,noRStudioGD = TRUE)
  #set_notebook_plot_size(20,20)
  #     png("gene_gene_interaction.png",width=3000,height=3000,res=250)
  logPInt <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {
    f<- try(fisher.test(genomicData[,i], genomicData[,j]), silent=TRUE)
    if(class(f)=="try-error"){
      0
    } else{
      ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))}}
  ))
  
  odds <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(genomicData[,i], genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  pairs <- sapply(1:ncol(genomicData), function(i) colMeans(genomicData * genomicData[,i], na.rm=TRUE))
  diag(logPInt) <- 0
  diag(odds) <- 1
  colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genomicData) 
  odds[odds<1e-3] = 1e-4
  odds[odds>1e3] = 1e4
  odds[10^-abs(logPInt) > 0.05] = 1
  logOdds=log10(odds)
  diag(logPInt) <- NA
  
  library(ComplexHeatmap)
  mat = logOdds
  rownames(mat) = columns
  colnames(mat) = columns
  #column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
  #row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))
  
  if (column_order){
    o = 1:length(columns) 
  } else{
    o = order(genomicGroups,-colSums(genomicData, na.rm=TRUE))#order(cmdscale(d, k=1))#h$order #c(h$order,(length(h$o rder) +1):ncol(interactions))
  }
  
  ## P values 
  P <- 10^-abs(logPInt[o,o])
  P_BH =  matrix(0, ncol = length(columns),nrow=length(columns)) 
  colnames(P_BH)<-columns 
  rownames(P_BH)<-columns
  P_FWER =  matrix(0, ncol = length(columns),nrow=length(columns)) 
  colnames(P_FWER)<-columns 
  rownames(P_FWER)<-columns
  
  P_BH[] = p.adjust(P, method="BH")
  P_FWER[] = p.adjust(P)
  P_BH[is.na(P_BH)] <- 999
  P_FWER[is.na(P_FWER)] <- 999
  
  P_sig =  matrix(0, ncol = length(columns),nrow=length(columns)) 
  colnames(P_sig)<-columns 
  rownames(P_sig)<-columns
  
  library(circlize)
  library(scales)
  
  
  col_fun = colorRamp2(c(-4, -3, -2, -1, 0-.Machine$double.eps, 0, 0+.Machine$double.eps,1, 2, 3), c("#8c510a","#bf812d", "#dfc27d", "#f6e8c3", "#f6e8c3","white", "#c7eae5","#c7eae5", "#80cdc1", "#35978f")) 
  p = Heatmap(mat, name = "Odds Ratio", cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = c("left"), column_names_side = c("top"), 
              heatmap_legend_param = list(
                at = c(-4, -3, -2, -1, 0, 1, 2),
                labels = c("0.0001", "0.001", "0.01", "0.1", "0", "10", "100"),
                title = "Odds Ratio",
                legend_height = unit(4, "cm"),
                title_position = "topleft"
              ),
              col = col_fun, rect_gp = gpar(col = "white", lwd = 2),
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(P_FWER[i, j] < 0.05) {
                  grid.text("\U2217", x, y)
                } else if(P_BH[i, j] < 0.1) {
                  grid.text("\U00B7", x, y)
                }
              })
  
  png("co-mut_hotspot_all_genes.png",width=13.5,height=12,units="in",res=1200)
  print(p)
  dev.off()
  pdf("co-mut_hotspot_all_genes.pdf",width=13.5,height=12)
  print(p)
  dev.off()
  print(p)
  
}


## load mutations files from Tazi et al, Papaemmanuil et al 
df_tazi = read.table("tazi_aml_prognosis_updated.tsv")
df_papaemme_nejm = read.table("df_nejm_personalization.tsv")
df_papaemme_nejm_add = df_papaemme_nejm %>% select("CEBPA_mono","CEBPA_bi","add_8","complex", "FLT3_other" ,"FLT3_TKD", "ITD","SRSF2","NRAS","add_11","add_13","add_21","add_22","del_12","del_17","del_18","del_20","del_3","del_4","del_5","del_7","del_9","inv_16","inv_3","minusy","t_15_17","t_6_9","t_8_21","t_9_11","t_9_22","t_v_11"    ) %>% 
  tibble::rownames_to_column("sample_pd")

# df_tazi = df_tazi %>% dplyr::mutate(NRAS = case_when(NRAS_p.G12_13 + NRAS_p.Q61_62 > 0 ~ 1, T ~ 0))

## papaemmanuil NEJM mutations 
x <- getURL("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt")
papaemme_vcf <- read.table(text = x, sep = "\t", header = TRUE)

## Split by IDH hotspot 
papaemme_vcf  = papaemme_vcf %>% mutate(gene_specific = case_when(GENE == "IDH2" & AA_CHANGE %in% c( "p.R140Q", "p.R140G", "p.R140L", "p.R140W") ~ "IDH2_p.R140",
                                                                  GENE == "IDH2" & AA_CHANGE %in% c( "p.R172K", "p.R172W") ~ "IDH2_p.R172",
                                                                  GENE == "NRAS" & AA_CHANGE %in% c("p.G12A", "p.G12C", "p.G12D", "p.G12R", "p.G12S", "p.G12V", "p.G13C", "p.G13D", "p.G13R", "p.G13V") ~ "NRAS_p.G12_13"  ,
                                                                  GENE == "NRAS" & AA_CHANGE %in% c("p.Q61H", "p.Q61K", "p.Q61L", "p.Q61P", "p.Q61R") ~ "NRAS_p.Q61_62"  ,
                                                                  T ~ GENE )) %>% mutate(sample_pd = SAMPLE_NAME) %>% filter(sample_pd %in% df_papaemme_nejm_add$sample_pd)
df_papaemme = papaemme_vcf %>% select(sample_pd, gene_specific) %>% 
  mutate(value = 1) %>%  distinct(sample_pd, gene_specific, .keep_all = TRUE) %>% 
  pivot_wider(names_from = "gene_specific", values_from = value, values_fill = 0)

df_papaemme = df_papaemme %>% full_join(df_papaemme_nejm_add, by = "sample_pd")
## length(unique(df_papaemme$sample_pd)) = 1540 

#df_tazi = df_tazi %>% tibble::rownames_to_column("sample_pd")
merged_cohorts = plyr::rbind.fill(df_tazi, df_papaemme)
## length(unique(merged_cohorts$sample_pd)) = 3653

merged_cohorts_IDH = merged_cohorts %>% filter(IDH2_p.R140 + IDH2_p.R172 + IDH1 >=1)


# check extra columns 
colnames(merged_cohorts)[!(colnames(merged_cohorts) %in% colnames(df_tazi))]
colnames(merged_cohorts)[!(colnames(merged_cohorts) %in% colnames(df_papaemme))]
colnames(df_papaemme_nejm)[!(colnames(df_papaemme_nejm) %in% colnames(df_papaemme))]


tmp <- merged_cohorts
#tmp = tmp %>% mutate(FLT3_ITD = ITD)

library(emojifont)


## merged NRAS hotspots and add CEBPA 
tmp1 = tmp %>% mutate(NRAS = case_when(NRAS == 1 | NRAS_p.G12_13 == 1 | NRAS_p.Q61_62 == 1 ~ 1, T ~ 0)) %>% 
  mutate(CEBPA = case_when(CEBPA == 1 | CEBPA_mono == 1 | CEBPA_bi == 1 ~ 1, T ~ 0) ) %>% 
  mutate(U2AF1 = case_when(U2AF1 == 1 | U2AF1_p.S34 == 1 | U2AF1_p.Q157 == 1 ~ 1, T ~ 0) ) %>% 
  mutate(SRSF2 = case_when(SRSF2 == 1 | SFRS2 == 1 ~ 1,  T ~ 0) )  %>% 
  mutate(FLT3 = case_when(FLT3 == 1 | FLT3_TKD == 1 | ITD == 1 | FLT3_other == 1 ~ 1 , T ~ 0) ) %>%
  mutate(KMT2C = case_when(MLL3 ==1 | KMT2C == 1 ~ 1, T ~ 0)) %>% 
  mutate(KMT2D = case_when(MLL2 ==1 | KMT2D == 1 ~ 1, T ~ 0)) %>% 
  mutate(KMT2E = case_when(MLL5 ==1 | KMT2E == 1 ~ 1, T ~ 0)) %>% 
  select(-NRAS_p.G12_13, -NRAS_p.Q61_62, -CEBPA_mono, -CEBPA_bi, -U2AF1_p.S34, -U2AF1_p.Q157,-SFRS2 , -FLT3_TKD, -ITD, -FLT3_other, -MLL3, -MLL3, -MLL5)

to.remove <- c(colnames(tmp1)[!(colnames(tmp1) %in% colnames(df_tazi))],
               colnames(tmp1)[!(colnames(tmp1) %in% colnames(df_papaemme))] )
to.remove = to.remove[! to.remove %in% c('FLT3', 'CEBPA', 'U2AF1', 'SRSF2', "NRAS", "KMT2C", "KMT2D", "KMT2E")]

`%ni%` <- Negate(`%in%`)

tmp2 = subset(tmp1,select = names(tmp1) %ni% to.remove)



colnames(tmp1 )
# set order 
all_columns = colnames(tmp2 )[! colnames(tmp2 ) == "sample_pd"]
ordered_cols = c("IDH1",  "IDH2_p.R140",  "IDH2_p.R172",  "TET2", "DNMT3A","ASXL1", "NPM1", "SRSF2",  "U2AF1", "STAG2", "CEBPA", "NRAS", "KRAS" ,"FLT3", "KIT", "BCOR", "add_8", "complex", "TP53",  "inv_3" ,"inv_16", "t_15_17", "t_v_11","t_6_9","t_8_21", "t_9_11","t_9_22"   )
other_cols = all_columns[! all_columns %in% ordered_cols]
columns_order = c(ordered_cols, other_cols)

gene_gene_interactions_OR(tmp2,columns=columns_order,
                          genomic_groups="Genetics",column_order=T)
