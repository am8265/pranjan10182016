#!/usr/bin/env Rscript
#####Installing the DESeq package########(Uncomment the following 2 lines if you would like to Install DESEQ package)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq")
library("DESeq")
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
cds<-estimateDispersions(cds)
str(fitInfo(cds))
res<-nbinomTest(cds, pair1, pair2)
###Load the TSV file containing annotations######
ptrich.allgenes = read.table(file = gene.annotation, header = TRUE, row.names = 1)
##select for genes at padjusted value of < 10% FDR####
res.sig<-res[ res$padj < 0.1, ]
if (is.na(res.sig[2:nrow(res.sig), ])) {#checking if ALL elements of the data frame is NA or not 
  print("res.sig is NA")
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
    res.sig.upreg = res.sig.upreg[order(- res.sig.upreg$foldChange, - res.sig.upreg$baseMean), ]
    res.sig.upreg100 = res.sig.upreg[c(1:100), ]#Top 100 upregulated genes
    ##remove all NA values!##############
    res.sig.upreg100=res.sig.upreg100[complete.cases(res.sig.upreg100), ]
    upreg100 = rownames(res.sig.upreg100)
    des = as.character(ptrich.allgenes[res.sig.upreg100$id, c(2)])
    res.sig.upreg100$description = des
    fname = paste(pair1, pair2,"sig", "upreg100","tsv",sep = ".")
    write.table(res.sig.upreg100, file = fname)
    sub_command = paste("sed -i.bak 's/_/ /g'",fname, sep = ' ')
    system(sub_command)
    #####variance stabilised data########
    cds.blind = estimateDispersions(cds,method="blind")
    vsd.blind = varianceStabilizingTransformation(cds.blind)
    fname.vsd.up.hmap = paste(pair1, pair2,"up.vsd", "hmap", "jpg", sep = ".")
    fname.norm.up.hmap=paste(pair1, pair2,"up.norm", "hmap", "jpg", sep = ".")
    fname.vsd.dwn.hmap = paste(pair1, pair2,"dwn.vsd", "hmap", "jpg", sep = ".")
    fname.norm.dwn.hmap=paste(pair1, pair2,"dwn.norm", "hmap", "jpg", sep = ".")
    library("RColorBrewer")
    library("gplots")
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    if (nrow(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]))>=2 & ncol(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]))>=2) {
      ######HeatMap for Upregulated VSD genes#####
      jpeg(fname.vsd.up.hmap)
      heatmap.2(as.matrix(exprs(vsd.blind)[as.numeric(upreg100), ]),col = hmcol, trace="none", margin=c(10, 6))
      dev.off()
      #######Heatmap for upregulated NORMALIZED genes################
      jpeg(fname.norm.up.hmap)
      heatmap.2(counts(cds,normalized=TRUE)[as.numeric(upreg100), ],col = hmcol, trace="none", margin=c(10, 6))
      dev.off()
    }  else {
         print("Due to Dimensionality problems we cannot produce a HeatMap!")
    }
    if (nrow(as.matrix(exprs(vsd.blind)[as.numeric(dwnreg100), ]))>=2 & ncol(as.matrix(exprs(vsd.blind)[as.numeric(dwnreg100), ]))>=2) {
      ######HeatMap for Downregulated VSD genes########
      jpeg(fname.vsd.dwn.hmap)
      heatmap.2(exprs(vsd.blind)[as.numeric(dwnreg100),],col = hmcol, trace="none", margin=c(10, 6))
      dev.off()
    ########Heatmap for downregulated Normalized genes############
      jpeg(fname.norm.dwn.hmap)
      heatmap.2(counts(cds,normalized=TRUE)[as.numeric(dwnreg100),],col = hmcol, trace="none", margin=c(10, 6))
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
