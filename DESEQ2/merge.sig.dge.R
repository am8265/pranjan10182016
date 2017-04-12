!#/usr/bin/env Rscript
#######The script reads "list.txt"-file containing relative paths to .csv files of DGE output of DESEQ2 run.
######The csv files of each CONDITION PAIR are merged to create 1)data.frames and .csv file
  #####containing all the signifcant (padj<0.05) DEGS across all 28 conditional pairs2)
  ###A subset of (1) containing DEGs which are ATLEAST 2 fold up/down in ATLEAST 1 conditional pair
setwd("/Users/ayanmalakar/pranjan10182016/DESEQ2")
filePath=getwd()
DGE.file.list<-read.table(file="list.txt", stringsAsFactors = F)#modify the list.txt as needed to provide
#28 conditional pairs
geneList=character()
geneList.core=character()
genes.of.int<-read.csv(file.path(filePath,DGE.file.list[,][1]), 
         header = T, row.names = 1)
geneList.1<-sort(row.names(genes.of.int))#gene list from the first .csv file of Gene of Interest
c=0
#parsing out the column names!
condition.pair<-strsplit(DGE.file.list[,][1], split=c("/"))[[1]][2]
condition.pair<-strsplit(condition.pair,split="csv")[[1]][1]
condition.pair<-paste(strsplit(condition.pair,split='\\.')[[1]],
                      collapse = '')
#####################
df<-data.frame(row.names = geneList.1)
df[condition.pair]<-genes.of.int[geneList.1,]$log2FoldChange####No threshold...
for (dge in DGE.file.list[2:28,1]) {
  
  c=c+1
  genes.of.int<-read.csv(file.path(filePath,dge), header = T, row.names = 1)
  geneList<-sort(row.names(genes.of.int))
  #geneList.union<-union(geneList,)
  #parsing out the column names!
  condition.pair<-strsplit(dge, split=c("/"))[[1]][2]
  condition.pair<-strsplit(condition.pair,split="csv")[[1]][1]
  condition.pair<-paste(strsplit(condition.pair,split='\\.')[[1]],
                        collapse = '')
  ##############################
  df1<-data.frame(row.names = geneList)
  df1[condition.pair]<-genes.of.int[geneList,]$log2FoldChange
  df<-merge(df,df1,by = "row.names", all=T)
  rownames(df)<-df[,1]
  df<-df[,-1]
  ###We comput geneList.core for Ctl vs. ALL conditions...
  # if (c==1) { 
  #   geneList.core<-intersect(geneList.1, row.names(genes.of.int))
  # } else {
  #   geneList.core<-intersect(geneList.core,row.names(genes.of.int))
  # }
}

ptrich.allgenes = read.table(file = "Ptrichocarpa_210_v3.0.allgenes.mod1.txt", header = TRUE, row.names = 1)
des<-as.character(ptrich.allgenes[rownames(df),]$description)
df["description"]<-des

fname = "condition.pairs.sig.dge.csv"
df.final<-cbind(gid=as.character(rownames(df)),df)
rownames(df.final)<-NULL
saveRDS(df.final, file="df.final.rds")
write.csv(df.final, file = fname)
sub_command = paste("sed -i.bak 's/_/ /g'",fname, sep = ' ')
system(sub_command)

########subset the  df.final for atleast 1 conditionpair -1>log2foldchange>=1######
df.final.mod<-df.final
rownames(df.final.mod)<-df.final$gid
head(df.final.mod)
df.final.mod$gid<-NULL
fold2<-which((df.final.mod>=1|df.final.mod<=-1),
             arr.ind = T)###select the genes with atleast 1 2foldchangedf.
df.2fold<-df.final.mod[rownames(fold2),]
fname = "condition.pairs.sig.2fold.dge.csv"
df.2fold<-cbind(gid=as.character(rownames(df.2fold)),df.2fold)
rownames(df.2fold)<-NULL
saveRDS(df.2fold, file="df.2fold.rds")
write.csv(file = fname, x = df.2fold)
#write.csv(df.final, file = fname)
sub_command = paste("sed -i.bak 's/_/ /g'",fname, sep = ' ')
system(sub_command)
# geneList.union<-union(geneList)
#########Creating the tabular conditional pairs file#########
