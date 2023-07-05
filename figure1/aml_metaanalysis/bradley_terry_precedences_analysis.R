options(repr.plot.width=10, repr.plot.height=6)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
source("my_bradley_terry.R")
library(RCurl)

# Load Tazi et al data
x <- getURL("https://raw.githubusercontent.com/papaemmelab/Tazi_NatureC_AML/main/data/genetic_files_main.tsv")
tazi_vcf <- read.csv(text = x)

# Load Papaemmanuil et al data 
x <- getURL("https://raw.githubusercontent.com/gerstung-lab/AML-multistage/master/data/AMLSG_Genetic.txt")
papaemme_vcf <- read.table(text = x, sep = "\t", header = TRUE)
papaemme_vcf = papaemme_vcf %>% mutate(sample_pd = SAMPLE_NAME, tum_vaf = VAF, tum_depth = TUM_DEPTH, chr = CHR, protein = AA_CHANGE, gene = GENE)

# Add patient sex to mutation file
clinical_papaemme <- read.table("full_data_validation.tsv", header=T, sep=" ", stringsAsFactors=F) # clinical papaemme
clinical_papaemme$sample_pd <- clinical_papaemme$Sample_ID
clinical_papaemme_subset = clinical_papaemme %>% select(sample_pd, Gender)
papaemme_vcf = papaemme_vcf %>% left_join(clinical_papaemme_subset) # Gender 0:723 female // 1:817 male 
papaemme_vcf$gene[papaemme_vcf$gene == 'SFRS2'] <- 'SRSF2'

clinical_tazi = read.table("tazi_aml_prognosis_updated.tsv")
clinical_tazi$sample_pd = rownames(clinical_tazi)
clinical_tazi$Gender = clinical_tazi$gender
clinical_tazi_subset = clinical_tazi %>% select(sample_pd, Gender)
tazi_vcf = tazi_vcf %>% left_join(clinical_tazi_subset) # Gender 0:929 female // 1:1184 male 

# Merge VCFs 
papaemme_vcf_tomerge = papaemme_vcf %>% select(sample_pd, Gender, chr, gene, protein, tum_vaf, tum_depth)
tazi_vcf_tomerge = tazi_vcf %>% select(sample_pd, Gender, chr, gene, protein, tum_vaf, tum_depth)
merged_vcf = rbind(tazi_vcf_tomerge, papaemme_vcf_tomerge)


# Inputs:
#   - maf: mutation file in maf format
#   !! in maf we expect columns: "TARGET_NAME","GENE","DEPTH","VAF","TCF"
#   - vec genes: vector of gene names for which precedences are to be evaluated
#   - ftest.threshold: fisher test p-value threshold for significant clonal heterogeneity in a sample
#   - num.prec.cutoff: minimum number of precedences to include gene in BT model
#   - reference gene for BT model --> is null then chosen from maf file



# Remove PTD and Missing VAF
merged_vcf <- merged_vcf %>% filter(!is.na(tum_vaf)) # %>% filter(type!="PTD") 

# Adjust VAF
merged_vcf$tum_vaf = as.numeric(merged_vcf$tum_vaf)
merged_vcf$tum_depth = as.numeric(merged_vcf$tum_depth)
merged_vcf <- merged_vcf %>% filter(!is.na(tum_vaf)) # %>% filter(type!="PTD") 
merged_vcf$tum_vaf = merged_vcf$tum_vaf/100

merged_vcf <- merged_vcf %>% mutate(vaf_corrected = case_when((chr == "X") & (Gender == 1) ~ tum_vaf/2, #  correct chr X for MALE 
                                                              (tum_vaf > 0.6) ~ tum_vaf/2, # assume LOH is VAF > 0.6 --> pretty rough but okay
                                                              TRUE ~ tum_vaf
))


# 
merged_maf <- merged_vcf %>% select(sample_pd, gene, tum_depth, tum_vaf, vaf_corrected)
#merged_maf$tum_vaf = merged_maf$tum_vaf/100
#merged_maf$vaf_corrected = merged_maf$vaf_corrected/100
merged_maf$tcf <- pmin(.99, 2*merged_maf$vaf_corrected)
colnames(merged_maf) = c("TARGET_NAME","GENE","DEPTH","VAF", "VAF_CORRECT", "TCF")

# ~~~~~~~~~
# For IDH2: 
patient.idh2 <- merged_maf %>% filter(GENE %in% c("IDH2"))  %>% pull(TARGET_NAME) %>% unique
merged_maf_IDH2_all =  merged_maf %>% filter(TARGET_NAME %in% patient.idh2)

t <- rev(sort(table(merged_maf_IDH2_all$GENE)))
genes_include <- names(t)[t>=15] # Select genes mutated in more than XX times in the IDH2 subset
genes_include = genes_include[genes_include != "TET2"] # exclude TET2 - mutually exclusive

num.prec.cutoff = 10
ref_gene = "IDH2"

c.res <- GetPrecedenceBradleyTerry(maf=merged_maf_IDH2_all, vec.genes=genes_include, 
                                   ftest.threshold=0.05, num.prec.cutoff=num.prec.cutoff, refBT=ref_gene)
c.plot <- PlotBT(resbt=c.res$resbt,maf=merged_maf_IDH2_all) #,colfill=dallcol[dallcol$comp==cc,"col"])
ggarrange(c.plot$ggvaf, c.plot$ggbt, widths = c(1,1.4))
ggsave("IDH2_BT_tazi_papaemme.png", width = 7, height = 8, useDingbats=FALSE)
ggsave("IDH2_BT_tazi_papaemme.pdf", width = 7, height = 8, useDingbats=FALSE)

# ~~~~~~~~~
# For IDH1: 

patient.idh1 <- merged_maf %>% filter(GENE %in% c("IDH1"))  %>% pull(TARGET_NAME) %>% unique
merged_maf_IDH1_all =  merged_maf %>% filter(TARGET_NAME %in% patient.idh1)

t <- rev(sort(table(merged_maf_IDH1_all$GENE)))
genes_include <- names(t)[t>=15] # Select genes mutated in more than XX times in the IDH2 subset

num.prec.cutoff = 10
ref_gene = "IDH1"

c.res <- GetPrecedenceBradleyTerry(maf=merged_maf_IDH1_all, vec.genes=genes_include, 
                                   ftest.threshold=0.05, num.prec.cutoff=num.prec.cutoff, refBT=ref_gene)
c.plot <- PlotBT(resbt=c.res$resbt,maf=merged_maf_IDH1_all) #,colfill=dallcol[dallcol$comp==cc,"col"])

pdf("IDH1_BT_text.pdf", width = 7, height = 8, useDingbats=FALSE)
ggarrange(c.plot$ggvaf, c.plot$ggbt, widths = c(1,1.4))

dev.off()

ggsave("IDH1_BT_tazi_papaemme.png", width = 7, height = 8)
ggsave("IDH1_BT_tazi_papaemme.pdf", width = 7, height = 8, device = cairo_pdf)



