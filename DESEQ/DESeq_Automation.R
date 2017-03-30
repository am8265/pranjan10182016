library("DESeq")
for (item in list.files(pattern = "*.csv")){
  sep_names=strsplit(item,split="\\.")[[1]][c(1,2)]
  dir_name=paste(sep_names,sep="",collapse = ".")
  print(dir_name)
  dir.create(dir_name)
  sample_pair=read.csv(item,row.names=1)
  sample_pair=sample_pair[rowSums(sample_pair)>0,]
  pair1=sep_names[c(1)]
  pair2=sep_names[c(2)]
  if (substr(pair1,3,4)=='-') {
    pair1_occur=grepl(substr(pair1,1,2),colnames(sample_pair),perl = "TRUE")
  } else {
      pair1_occur=grepl(pair1,colnames(sample_pair))
  }
  pair1_no=length(pair1_occur[pair1_occur==TRUE])
  pair2_no=ncol(sample_pair)-pair1_no
  condition=c(rep.int(pair1,pair1_no),rep.int(pair2,pair2_no))
  cds = newCountDataSet(sample_pair,condition)
  cds = estimateSizeFactors( cds )
  sizeFactors( cds )
  path=paste(getwd(),dir_name,sep="/")
  fname_norm=paste(pair1,pair2,"norm","csv",sep=".")
  norm_file=counts( cds, normalized=TRUE )
  write.csv(norm_file,paste(path,fname_norm,sep='/'))
  cds = estimateDispersions( cds)
  str(fitInfo(cds))
  path=paste(getwd(),dir_name,sep="/")
  fname_disp=paste(pair1,pair2,"displot","jpg",sep=".")
  path_fname_disp=paste(path,fname_disp,sep='/')
  jpeg(path_fname_disp)
  plotDispEsts( cds )
  dev.off()
  head( fData(cds) )
  res = nbinomTest(cds,pair1,pair2)
  fname_MA=paste(pair1,pair2,"MA","jpg",sep=".")
  path_fname_MA=paste(path,fname_MA,sep = '/')
  jpeg(path_fname_MA)
  plotMA(res)
  dev.off()
  fname_hist=(paste(pair1,pair2,"pval","jpg",sep = "."))
  path_fname_hist=paste(path,fname_hist,sep="/")
  jpeg(path_fname_hist)
  hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
  dev.off()
  resSig = res[ res$padj < 0.1, ]
  head( resSig[ order(resSig$pval), ] )
  resSig_order_pval=resSig[ order(resSig$pval), ]
  resSig_dwnreg= resSig[order( resSig$foldChange, -resSig$baseMean), ]
  resSig_upreg= resSig[ order( -resSig$foldChange, -resSig$baseMean ), ]
  fname_dwnreg=paste(pair1,pair2,"res.dwnreg","csv",sep=".")
  fname_upreg=paste(pair1,pair2,"res.upreg","csv",sep=".")
  write.csv(resSig_dwnreg,paste(path,fname_dwnreg,sep='/'))
  write.csv(resSig_upreg,paste(path,fname_upreg,sep='/'))
}
  
  
  
}
  #sep_names=strsplit(item,split="\\.")[[1]][c(1,2)]
#dir.create(dir_name)