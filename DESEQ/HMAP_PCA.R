#!/usr/bin/env Rscript
#####Installing the DESeq package########(Uncomment the following 2 lines if you would like to Install DESEQ package)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
library("DESeq")
setwd("/Users/ayanmalakar/pranjan10182016/DESEQ")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least 2 arguments must be supplied (input file).n", call. = FALSE)
}
item = args[1]
gene.annotation = args[2] 
sep.names<-strsplit(item,split="\\.")[[1]][c(1, 2)]#split a string into a character vector using the pattern
##works like a split function in python
sample.pair<-read.csv(item, row.names = 1)
sample.pair<-sample.pair[rowSums(sample.pair) > 0, ]#Collect genes with counts >0
pair1<-sep.names[c(1)]
pair2<-sep.names[c(2)]
#####loading the cumulative pvalues if they exist#####
if (file.exists("pval.FC0.cumulative.rds")) {
  pval.FC0.cumulative=readRDS("pval.FC0.cumulative.rds")
} else {
  pval.FC0.cumulative=vector("numeric")
}

if (file.exists("padj.FC0.cumulative.rds")){
  padj.FC0.cumulative=readRDS("padj.FC0.cumulative.rds")
} else {
  padj.FC0.cumulative=vector("numeric")
}
if (substr(pair1, 3, 4) == '-')  {#inspect the 3rd character
  pair1.occur<-grepl(substr(pair1, 1, 2), colnames(sample.pair), perl = "TRUE")
} else  {
    pair1.occur<-grepl(pair1,colnames(sample.pair)) 
}
pair1.no<-length(pair1.occur[pair1.occur == TRUE])
pair2.no<-ncol(sample.pair) - pair1.no
condition<-factor(c(rep.int(pair1, pair1.no),rep.int(pair2, pair2.no)))
cds<-newCountDataSet(sample.pair, condition)
cds<-estimateSizeFactors(cds)
sizeFactors(cds)
write.csv(as.data.frame(counts(cds,normalized = TRUE)),file=paste(pair1,pair2,"normalised","csv",sep="."))
cds<-estimateDispersions(cds)
str(fitInfo(cds))
res<-nbinomTest(cds, pair1, pair2)
####cumulative pval and padj#######
jpeg(paste(pair1,pair2,"his.pval.FC0.jpeg",sep = "."))
pval.FC0=ifelse(res$foldChange==0,res$pval,NA)
pval.FC0=pval.FC0[!is.na(pval.FC0)]
saveRDS(pval.FC0,paste(pair1,pair2,"pval.FC0.rds",sep="."))
pval.FC0.cumulative=c(pval.FC0.cumulative,pval.FC0)
saveRDS(pval.FC0.cumulative,"pval.FC0.cumulative.rds")
if (length(pval.FC0)!=0) {
  hist(pval.FC0, breaks = 100, col = "skyblue", border = "slateblue", main=" ")
}
dev.off()
jpeg(paste(pair1,pair2,"his.padj.FC0.jpeg",sep = "."))
padj.FC0=ifelse(res$foldChange==0,res$padj,NA)
padj.FC0=padj.FC0[!is.na(padj.FC0)]
saveRDS(padj.FC0,paste(pair1,pair2,"padj.FC0.rds",sep="."))
padj.FC0.cumulative=c(padj.FC0.cumulative,padj.FC0)
saveRDS(padj.FC0.cumulative,"padj.FC0.cumulative.rds")
if (length(padj.FC0)!=0) {
  hist(padj.FC0, breaks=100, col = "skyblue", border = "slateblue", main=" ")
}
dev.off()
saveRDS(res, paste(pair1,pair2,"res","rds",sep = "."))###save the RDS file here 
###Load the TSV file containing annotations######
ptrich.allgenes = read.table(file = gene.annotation, header = TRUE, row.names = 1)
##select for genes at padjusted value of < 10% FDR####
res.sig<-res[ res$padj < 0.1, ]
if (is.na(res.sig[2:nrow(res.sig), ])) {#checking if ALL elements of the data frame is NA or not 
  print("res.sig is NA")
  saveRDS(res.sig, paste(pair1,pair2,"rds",sep="."))
  
}  else {
    res.sig.dwnreg = res.sig[res.sig$foldChange<1,]#Select only genes with fold change <1
    #among the genes with fold change<1,order them from lowest value to highest (break a tie with higher gene expression across samples)
    res.sig.dwnreg = res.sig.dwnreg[order(res.sig.dwnreg$foldChange, -res.sig.dwnreg$baseMean), ]
    res.sig.dwnreg100 = res.sig.dwnreg[c(1:100), ]##Top 100 significantly Downregulated genes
    #remove all NA values############## 
    res.sig.dwnreg100 = res.sig.dwnreg100[complete.cases(res.sig.dwnreg100), ]
    dwnreg100 = rownames(res.sig.dwnreg100)
    des = as.character(ptrich.allgenes[res.sig.dwnreg100$id, c(2)])
    res.sig.dwnreg100$description = des
    fname = paste(pair1, pair2,"sig", "dwnreg100","tsv", sep = ".")
    write.table(res.sig.dwnreg100, file = fname)
    sub_command = paste("sed -i.bak 's/_/ /g'",fname, sep = ' ')
    system(sub_command)
    ####Significantly Upregulated genes##########
    res.sig.upreg = res.sig[res.sig$foldChange>1,]
    res.sig.upreg = res.sig.upreg[complete.cases(res.sig.upreg), ]
    res.sig.upreg = res.sig.upreg[order(- res.sig.upreg$foldChange, - res.sig.upreg$baseMean), ]
    res.sig.upreg100 = res.sig.upreg[c(1:100), ]#Top 100 upregulated genes
    ##Remove all NA values among significantly upregulated genes!##############
    res.sig.upreg100=res.sig.upreg100[complete.cases(res.sig.upreg100), ]
    upreg100 = rownames(res.sig.upreg100)
    des = as.character(ptrich.allgenes[res.sig.upreg100$id, c(2)])
    res.sig.upreg100$description = des
    fname = paste(pair1, pair2,"sig", "upreg100","tsv",sep = ".")
    write.table(res.sig.upreg100, file = fname)
    sub_command = paste("sed -i.bak 's/_/ /g'",fname, sep = ' ')
    system(sub_command)
    #######Combining the significantly expressed genes and saving it into R data file#
    res.sig.combine = res.sig[res.sig$foldChange>1 | res.sig$foldChange<1,]#we are skipping where fold change =1
    write.csv(res.sig.combine, file = paste(pair1, pair2,"sig","csv",sep = "."))
    rownames(res.sig.combine)=NULL
    saveRDS(res.sig.combine, paste(pair1,pair2,"rds",sep="."))
    ######Combining the 2fold significantly expressed genes and saving it into R data file######
    res.sig.combine.2fold = res.sig[res.sig$log2FoldChange>=1 | res.sig.combine$log2FoldChange<=-1, ]
    write.csv(res.sig.combine.2fold, file = paste(pair1, pair2,"sig.2foldChange","csv", sep="."))
    rownames(res.sig.combine.2fold)=NULL
    saveRDS(res.sig.combine.2fold, paste(pair1,pair2,"2fold","rds",sep="."))
    #####variance stabilised data########
    cds.blind = estimateDispersions(cds,method="blind")
    vsd.blind = varianceStabilizingTransformation(cds.blind)
    write.csv(exprs(vsd.blind),file=paste(pair1,pair2,"variance.stabilised","csv",sep="."))
    fname.vsd.up.hmap = paste(pair1, pair2,"up.vsd", "hmap", "jpg", sep = ".")
    fname.norm.up.hmap=paste(pair1, pair2,"up.norm", "hmap", "jpg", sep = ".")
    fname.vsd.dwn.hmap = paste(pair1, pair2,"dwn.vsd", "hmap", "jpg", sep = ".")
    fname.norm.dwn.hmap=paste(pair1, pair2,"dwn.norm", "hmap", "jpg", sep = ".")
    library("RColorBrewer")
    library("gplots")
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    if (nrow(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]))>=2 & ncol(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]))>=2) {
      ######HeatMap for Upregulated VSD genes#####
      jpeg(fname.vsd.up.hmap, height=1000, width=1000)
      heatmap.2(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]),col = hmcol, trace="none", margin=c(6, 6),srtCol=15, offsetCol =-1.5,cexCol=1.1)
      dev.off()
      #######Heatmap for upregulated NORMALIZED genes################
      jpeg(fname.norm.up.hmap, height=1000, width=1000)
      heatmap.2(counts(cds,normalized=TRUE)[as.numeric(upreg100), ],col = hmcol, trace="none", margin=c(6, 6),srtCol=15, offsetCol =-1.5,cexCol=1.1)
      dev.off()
    }  else {
         print("Due to Dimensionality problems we cannot produce a HeatMap!")
    }
    if (nrow(as.matrix(exprs(vsd.blind)[as.numeric(dwnreg100), ]))>=2 & ncol(as.matrix(exprs(vsd.blind)[as.numeric(dwnreg100), ]))>=2) {
      ######HeatMap for Downregulated VSD genes########
      jpeg(fname.vsd.dwn.hmap,height=1000, width=1000)
      heatmap.2(exprs(vsd.blind)[as.numeric(dwnreg100),],col = hmcol, trace="none", margin=c(6, 6),srtCol=15, offsetCol =-1.5,cexCol=1.1)
      dev.off()
    ########Heatmap for downregulated Normalized genes############
      jpeg(fname.norm.dwn.hmap, height=1000, width=1000)
      heatmap.2(counts(cds,normalized=TRUE)[as.numeric(dwnreg100),],col = hmcol, trace="none", margin=c(6, 6),srtCol=15, offsetCol =-1.5,cexCol=1.1)
      dev.off()
    } else {
        print("Due to Dimensionality problems we cannot produce a HeatMap!")
    }
  ########PCA plots based on VSD##############
  fname.vsd.pca = paste(pair1, pair2, "pca","jpg",sep=".")#change this line of code
  jpeg(fname.vsd.pca)
  print(plotPCA(vsd.blind,intgroup=c("condition")))
  dev.off()
}
