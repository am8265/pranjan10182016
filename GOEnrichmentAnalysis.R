#!/usr/bin/env Rscript
######Run this script as ./GOEnrichmentAnalysis.R <file path to count matrix>
###settin up the working directory#####
setwd("/Users/ayanmalakar/pranjan10182016/GOEXPRESS/")
#Importing all the Annotation files#
ptrich.allgenes<-read.table(file='../Ptrichocarpa_210_v3.0.allgenes.mod1.txt',header = TRUE)
ptrich.allGO<-read.csv(file="../Ptrichocarpa_210_v3.0.allGO.txt",header = TRUE)
ptrich.GOgenes<-read.table(file='../Ptrichocarpa_210_v3.0.gi_go.txt',header = TRUE)
library(GOexpress)
#CO4.Ctl=read.csv(file='CO4.Ctl.csv', header = TRUE, row.names = 1)
library("Biobase")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least 1 arguments must be supplied (input file).n", call. = FALSE)
}
exprsFile = args[1]
exprs <- as.matrix(read.csv(exprsFile, header=TRUE,row.names=1, as.is = TRUE))
exprs<-exprs[rowSums(exprs)>0,]#deleting all 0 rows if such rows exist!
sep.names<-strsplit(exprsFile,split="\\.")[[1]][c(3, 4)]#parsing the conditions
pair1<-strsplit(sep.names[c(1)],split = "/")[[1]][2]
pair2<-sep.names[c(2)]
if (substr(pair1, 3, 4) == '-')  {#inspect the 3rd character
  pair1.occur<-grepl(substr(pair1, 1, 2), colnames(exprs), perl = "TRUE")
} else  {
  pair1.occur<-grepl(pair1,colnames(exprs)) 
}
pair1.no<-length(pair1.occur[pair1.occur == TRUE])
pair2.no<-ncol(exprs) - pair1.no
######Crearting the ExpressionSet Class########
pData<-data.frame(row.names = colnames(exprs),treatment = c(rep(pair1,pair1.no),rep(pair2,pair2.no)))
phenoData<-new("AnnotatedDataFrame",data=pData)
exprsSet<-ExpressionSet(assayData = exprs,phenoData = phenoData)
set.seed(4543)#set random seed for reproducibility
exprsSet.results<-GO_analyse(eSet = exprsSet, f="treatment", all_genes = ptrich.allgenes, all_GO = ptrich.allGO, GO_genes = ptrich.GOgenes)
saveRDS(exprsSet.results,file=paste(pair1,pair2,"results.rds",sep="."))
######Some Diagnostics###########
head(exprsSet.results$GO[,c(1:5,7)])#Ranked table of GO terms (subset)
head(exprsSet.results$genes[,c(1:3)])#Ranked table of genes(subset)
head(exprsSet.results$mapping)#gene to gene ontology mapping table
############Permutation-based P-value for ontologies##################
exprsSet.results.pVal=pValue_GO(result = exprsSet.results, N=100)#Increase to atleast N=1000 permutations while running on newton
#########Filtering exprsSet.results.pVal#########
BP.5<-subset_scores(result=exprsSet.results.pVal,namespace="biological_process",total=5, p.val=0.05)
MF.10<-subset_scores(result=exprsSet.results.pVal,namespace="molecular_function", total=10, p.val=0.05)
CC.15<-subset_scores(result=exprsSet.results.pVal,namespace="cellular_component",total=15,p.val=0.05)
#############Top Ranking GO Terms based on BP.5 , MF.10 and CC.15#######
head(BP.5$GO)
head(MF.10$GO)
head(CC.15$GO)
####################Hierarchial clustering of samples based on gene expression associated with a GO Terms######
######TOP GO Terms from BP.5#####
dir.create(paste(pair1,pair2,sep = '_'))#creates a directory based on individual sample names
ontology=c("BP.5","MF.10","CC.15")
for (go in ontology){
  for (term in as.character((head(get(go)$GO))$go_id)){
    jpeg(paste(pair1,pair2,go,term,"jpg",sep='.'), height=1000, width=1000)
    heatmap_GO(go_id = term, result =get(go), gene_names=FALSE, eSet = exprsSet, cexRow = 0.4, cexCol = 1,cex.main = 1,main.Lsplit = 30)
    dev.off()
    file.rename(paste(pair1,pair2,go,term,"jpg",sep='.'),paste(pair1,"_",pair2,"/",paste(pair1,pair2,go,term,"jpg",sep='.'),sep=''))
  }
}
#jpeg(paste(pair1,pair2,), height=1000, width=1000)