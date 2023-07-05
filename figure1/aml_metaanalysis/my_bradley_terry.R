
library("BradleyTerry2")
GetPrecedenceBradleyTerry <- function(maf, vec.genes=NULL, ftest.threshold=0.01, num.prec.cutoff=10, refBT="SF3B1") {
  
  # Inputs:
  #   - maf: mutation file in maf format
  #   !! in maf we expect columns: "TARGET_NAME","GENE","DEPTH","VAF","TCF"
  #   - vec genes: vector of gene names for which precedences are to be evaluated
  #   - ftest.threshold: fisher test p-value threshold for significant clonal heterogeneity in a sample
  #   - num.prec.cutoff: minimum number of precedences to include gene in BT model
  #   - reference gene for BT model --> is null then chosen from maf file
  
  # Outputs:
  # - precedence: matrix of precedences
  # - precedenceLong: df of precedences
  # - precedenceFull: matrix of precedences before "num.prec.cutoff" cut
  # - patientList: list of patient
  
  # Preparation
  if (is.null(vec.genes)) {
    vec.genes <- unique(maf$GENE)
  }
  maf <- maf[maf$GENE %in% vec.genes,]
  vec.samples <- unique(maf$TARGET_NAME)[sapply(unique(maf$TARGET_NAME),function(p) length(unique(maf[maf$TARGET_NAME==p,"GENE"]))>=2)] # >= 2 mutated genes of interest
  maf <- maf[maf$TARGET_NAME %in% vec.samples,]
  if (is.null(refBT)) {
    refBT <- names(rev(sort(table(maf$GENE))))[2]
    print(refBT)
  }
  
  # Precendence Matrix
  # Row Gene Before Col Gene
  precedence <- matrix(0, nrow=length(vec.genes), ncol = length(vec.genes) , dimnames=list(vec.genes,vec.genes))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  plist <- list()
  for (ss in unique(maf$TARGET_NAME)) {
    
    list.pair <- list()
    
    ssmaf = maf[maf$TARGET_NAME==ss,]
    ssmaf$clonality_label <- NA # think about it
    
    for (i in 1:nrow(ssmaf)) {
      for (j in 1:nrow(ssmaf)) {
        
        if (i!=j) {
          # ~~~~~~~~~~ #
          # ~~~~ Get Matrix of significance precendence ~~~~ #
          # Null Hypothesis:
          # the two mutations are present in the same fraction of cells
          # if reject the null --> the two mutations are likely present in different factions of cells
          tcf1 = ssmaf[i,"TCF"]
          depth1 = ssmaf[i,"DEPTH"]
          
          tcf2 = ssmaf[j,"TCF"]
          depth2 = ssmaf[j,"DEPTH"]
          
          lesions1 = ssmaf[i,"GENE"]
          lesions2 = ssmaf[j,"GENE"]
          
          # Contingency table of expected number of reads that explained the TCF proportions
          m <- round(matrix(c(
            tcf1*depth1,
            depth1-tcf1*depth1, 
            tcf2*depth2,
            depth2-tcf2*depth2),ncol=2))
          
          #f <- try(fisher.test(m, alternative="greater")$p.value< 0.01 , silent=TRUE)
          ftest <- ( fisher.test(m, alternative="greater")$p.value < ftest.threshold )
          
          if (ftest) { # significant clonal heterogeneity
            if (tcf1 + tcf2 >= 1) { # pigeon-hole --> nested clone [not parallel clones] 
              # Apply precendence
              precedence[lesions1,lesions2] <- precedence[lesions1,lesions2] + 1
              list.pair <- c(list.pair, list(c(lesions1,lesions2)))
            }
          }
          # ~~~~~~~~~~ #
        }
      }
    }
    plist[[ss]] <- list.pair
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  precedence0 = precedence
  
  # Remove genes with too little informative precedence
  diag(precedence) <- 0
  sump <- apply(precedence,1,sum)
  precedence <- precedence[which(sump>=num.prec.cutoff),which(sump>=num.prec.cutoff)]
  
  # Transform into long df
  precedence.sf <- countsToBinomial(precedence)
  names(precedence.sf)[1:2] <- c("gene1", "gene2")
  # Apply Bradley Terry
  bModel <- BTm(cbind(win1, win2), gene1, gene2, ~ gene,
                id = "gene", data = precedence.sf, refcat=refBT)
  resbt <- BTabilities(bModel)
  resbt <- data.frame(resbt)
  resbt$ability<- as.numeric(-resbt$ability)
  resbt <- resbt[order(resbt$ability,decreasing=T),]
  resbt$gene <- rownames(resbt)
  resbt$gene <- factor(resbt$gene, levels=resbt$gene)
  resbt$number <- sump[as.vector(resbt$gene)]
  
  return(list(precedence=precedence, precedenceLong=precedence.sf, precedenceFull=precedence0,
              resbt=resbt,patientlist=plist
  ))
  
}


#t <- table(sapply(plist, length)>0)
#pie(t, labels=paste(t, c("clonal/NA","polyclonal")), col=c("grey","orange"))

library(ggridges)
PlotBT <- function(resbt,maf,colfill="grey") {
  
  # BT plot
  ggbt <- ggplot(resbt, aes(y=gene, x=ability)) + 
    geom_segment( aes(y=gene, yend=gene, x=ability-1.96*s.e., xend=ability+1.96*s.e.),color='darkgrey') +
    geom_point(size=2) + 
    geom_text(aes(x=ability-1.96*s.e.-0.3,y=gene,label = gene)) + 
    geom_text(aes(x=ability+1.96*s.e.+0.3,y=gene,label = paste0("n=",number))) + 
    theme_classic() + 
    theme(axis.line = element_blank(),
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          text=element_text(size=15)
    ) + #noytitle +
    annotate("segment", x = min(resbt$ability-1.96*resbt$s.e.),xend = max(resbt$ability+1.96*resbt$s.e.),
             y = 0, yend = 0, colour = "black", size=1, alpha=1, arrow=arrow()) + 
    expand_limits(y=c(-0.5,+0.5), x=-0.6) + xlab("Relative time")
  
  # VAF density plot
  btmaf <- maf[maf$GENE%in%resbt$gene,]
  btmaf$GENE = factor(btmaf$GENE, levels=resbt$gene)
  
  ggvaf <- ggplot(btmaf, aes(y=GENE,x=TCF)) + geom_density_ridges(fill=colfill,alpha=.8) + theme_classic() +
    theme( axis.line.y = element_blank(),text=element_text(size=15),axis.ticks.y = element_blank()) + #noytitle + 
    xlab("Adjusted VAF")
  
  return(list(ggvaf=ggvaf,ggbt=ggbt))
}
