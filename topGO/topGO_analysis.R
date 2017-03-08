#!/usr/bin/env Rscript
source("http://bioconductor.org/biocLite.R")
biocLite("topGO")
library(topGO)
library(ALL)
biocLite("ALL")
biocLite("GO.db")
########Change the readRDS file at every stage##############
res=readRDS('../DESEQ/Ctl.Lb-.res.rds')
geneList=res$padj
names(geneList)=res$id
topDiffGenes<-function(allScore){
  return(allScore<=0.05)
}
x<-topDiffGenes(geneList)
sum(x)
go_category=c(5,10,15)
names(go_category)=c("BP","MF","CC")
for (category in names(go_category)) {
  GOdata<-new("topGOdata",description=paste("GO analysis of Ctl vs Lb- samples",category,sep = ' '), allGenes=geneList,geneSel=topDiffGenes, gene2GO=geneID2GO, annot=annFUN.gene2GO, ontology=category, nodeSize=go_category[category])
  saveRDS(GOdata,paste("GOdata",category,"rds",sep="."))
  description(GOdata)
  #classical enrichment analysis by testing the over-representation of GO terms within the
  #group of differentially expressed genes
  resultFisher<-runTest(GOdata,algorithm="classic",statistic = "fisher")
  saveRDS(resultFisher,paste("resultFisher",category,"rds",sep="."))
  ##the enrichment using the Kolmogorov-Smirnov test..algorithm classic method
  resultKS <- runTest(GOdata, algorithm = "classic",statistic = "ks")
  saveRDS(resultKS,paste("resultKS.classic",category,"rds",sep = "."))
  #####algorithm elim method
  resultKS.elim <- runTest(GOdata, algorithm = "elim",statistic = "ks")
  saveRDS(resultKS.elim,paste("resultKS.elim",category,"rds",sep = "."))
  allRes<-GenTable(GOdata,classicFisher=resultFisher, classicKS=resultKS,elimKS=resultKS.elim, orderBy="elimKS", ranksOf="classicFisher")
  saveRDS(allRes,paste("allRes",category,"rds",sep = "."))
  # pValue.classic<-score(resultKS)
  # pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
  # gstat <- termStat(GOdata, names(pValue.classic))
  # gSize <- gstat$Annotated / max(gstat$Annotated) * 4
  # gCol <- colM
  # plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize, col = gCol)
  png(paste("subgraph",category,"png",sep="."),width=400,height = 350, res)
  showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
  dev.off()
}


