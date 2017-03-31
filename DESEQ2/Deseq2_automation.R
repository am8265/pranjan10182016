#!/usr/bin/env Rscript
#####loading the dependencies########
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
filePath<-getwd()
file.list=list.files(pattern='*.csv')
for (item in file.list){
  countData<-as.matrix(read.csv(file.path(filePath,item),header = TRUE,row.names = 1))
  sep.names<-strsplit(item,split="\\.")[[1]][c(1, 2)]
  pair1<-sep.names[c(1)]
  pair2<-sep.names[c(2)]
  dir.create(file.path(filePath,paste(pair1,pair2,sep="_")),showWarnings = F)
  setwd(file.path(filePath,paste(pair1,pair2,sep="_")))
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
  print(plotPCA(rld, intgroup=c("condition")))
  dev.off()
  #######A More Customised PCA plot#####################
  
  ####Identifying possible outliers for improving the results and performing a PCA plot on them !!!!#####
  ####To be decided ##############################
  
  #######Performing Some Initial Sanity Checks########
  GeneCounts <- counts(dds)
  idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
  ##In a typical RNA-Seq experiment, there will be at least several thousand genes that are expressed in all samples. If the number of non–zero genes is very low, there is usually something wrong with at least some of the samples (e.g. the sample preparation was not done properly, or there was some problem with the sequencing).
  nonzero.genes=sum(idx.nz)/(dim(dds)[1])*100
  nonzero.genes
  atleast.one.zero.genes=100-nonzero.genes
  atleast.one.zero.genes
  
  ###Normalising the dataset############
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  ###Assessing the normalisation
  png(paste(pair1,pair2,"multidensity.png",sep="."), height= 5*300, width = 5*300, res=300, pointsize = 8)
  multidensity(counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  png(paste(pair1,pair2,"multiecdf.png",sep="."), height= 5*300, width = 5*300, res=300, pointsize = 8)
  multiecdf( counts(dds, normalized = T)[idx.nz ,],xlab="mean counts", xlim=c(0, 1000))
  dev.off()
  #####Plotting pairwise MeanAverage Plots#########
  # png(paste(pair1,pair2,"pairwiseMAs.png", sep="."), height= 5*300, width = 5*300, res=300, pointsize = 8)
  # MA.idx = t(combn(1:ncol(countData), 2))
  # for( i in  seq_along( MA.idx[,1])){ 
  #   MDPlot(counts(dds, normalized = T)[idx.nz ,], 
  #          c(MA.idx[i,1],MA.idx[i,2]), 
  #          main = paste( colnames(dds)[MA.idx[i,1]], " vs ",
  #                        colnames(dds)[MA.idx[i,2]] ), ylim = c(-3,3))
  # }
  # dev.off()
  
  ####Assessing for DGE AS IT IS i.e we are NOT REMOVING ANY OUTLIERS########
  ##Estimating Dispersions#####################
  dds<-estimateDispersions(dds)
  png(paste(pair1,pair2,"dispPlot.png",sep="."), height= 5*300, width = 5*300, res=300, pointsize = 8)
  plotDispEsts(dds)
  dev.off()
  #####Testing for DGE using nbinomWaldTest for less than 5% FDR...
  #####...independent filtering is automatically performed 
  dds<-nbinomWaldTest(dds)
  if (pair1=="Ctl"){
    dds_res<-results(dds, pAdjustMethod = "BH", contrast = c("condition",pair2, pair1), alpha = 0.05)
  }
  dds_res<-results(dds, pAdjustMethod = "BH", alpha=0.05)
  summary(dds_res)
  sum(dds_res$padj<0.05, na.rm = TRUE)#No. of genes with significant DE
  sum(dds_res$padj<0.05 & (dds_res$log2FoldChange>=1 |dds_res$log2FoldChange<=-1), na.rm = T)#No. of genes with significant DE and 
  ###fold change of >=2 or <=0.5!
  saveRDS(dds_res,paste(pair1,pair2,"dds_res.rds",sep = "."))
  png(paste(pair1,pair2,"MAplot.png",sep="."),height= 5*300, width = 5*300, res=300, pointsize = 8)
  plotMA(dds_res, alpha=0.05, main=ifelse(pair1!="Ctl",paste(pair1,"vs.",pair2,sep=""),paste(pair2,"vs",pair1)))
  dev.off()
  dds_res<-dds_res[order(dds_res$padj),]
  resSig05.2fold<-subset(dds_res, padj< 0.05 & log2FoldChange>=1|log2FoldChange<=-1)
  resSig05<-subset(dds_res, padj< 0.05 & log2FoldChange!=0)
  write.csv(as.data.frame(resSig05.2fold), file=paste(pair1,pair2,"2fold","csv",sep="."))
  write.csv(as.data.frame(resSig05), file=paste(pair1,pair2,"csv", sep = "."))
  #######Testing for DGE using Log Likelihood Ratio Test#######
  
  
  
  
  
  ####REANALYSING EVERYTHING AFTER REMOVING OUTLIERS!!!!!!###########
  #######How to remove and select for outliers ##########
  
}








