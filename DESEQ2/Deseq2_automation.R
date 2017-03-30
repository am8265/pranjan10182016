#####Dowloading the dependencies########
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Mm.eg.db)
setwd("/Users/ayanmalakar/pranjan10182016/DESEQ2")
file.list=list.files(pattern='*.csv')
for (item in file.list){
  countData<-as.matrix(read.csv(item,header = TRUE,row.names = 1))
  sep.names<-strsplit(item,split="\\.")[[1]][c(1, 2)]
  pair1<-sep.names[c(1)]
  pair2<-sep.names[c(2)]
  if (substr(pair1, 3, 4) == '-')  {#inspect the 3rd character
    pair1.occur<-grepl(substr(pair1, 1, 2), colnames(countData), perl = "TRUE")
  } else  {
    pair1.occur<-grepl(pair1,colnames(countData)) 
  }
  pair1.no<-length(pair1.occur[pair1.occur == TRUE])
  pair2.no<-ncol(countData) - pair1.no
  colData<-data.frame(colnames(countData), condition=c(rep(pair1,pair1.no),rep(pair2,pair2.no)))
  dds<-DESeqDataSetFromMatrix(countData = countData,colData = colData, design= ~ condition)
  dds<-dds[rowSums(counts(dds))>0,]#pre-filtering for genes with atleast 1 non-zero count
  dds<-DESeqDataSetFromMatrix(countData = countData,colData = colData, design= ~ condition)
  ###Sample-Sample Clustering HeatMaps and PCA plot
  rld <- rlogTransformation(dds, blind = TRUE)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  rownames(mat)<-colData(rld)$condition
  ###sample-sample clustering and visualising on Heat Map
  hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(255)
  png(paste(pair1,pair2,"sampleCluster","png",sep="."), height= 5*300, width = 5*300, res=300, pointsize = 8)
  heatmap.2(mat, trace="none", col = rev(hmcol),cexCol = 1.0,srtCol = 15)
  dev.off()
  #####PCAplot for identifying outliers#####
  png(paste(pair1,pair2,"pca","top500.png",sep="."),height= 5*300, width = 5*300, res=300, pointsize = 8)
  plotPCA(rld, intgroup=c("condition"))
  dev.off()
  #######A More customised PCA plot#####################
  
  ####Identifying possible outliers for improving the results#####
  ####To be decided ##############################
  
  ####Assessing for DGE AS IT IS########
  ####
  
  
  
  
  
  ####Assessing for DGE after removing outliers###########
  









