How to Interpret the files provided in each directory?

PS:The order of conditions/sample group are same as provided in the <pair1.pair2>.<extension> prefix for each files
 
<pair1.pair2>.csv-Count values for 2 conditions as obtained after running HTSEQ on tophat output and using custom scripts
<pair1.pair2>.csv was used as an input to DESeq Bioconductor package to analyse for Differential Gene Expression.
eg. CO4.Ctl.csv. Here pair1 means CO4 and pair2 means Ctl

<pair1.pair2>.DGE.csv-Significant Differentially Expressed Genes across 2 Conditions.FDR or padj.<0.1 was used for Multiple testing correction to find significant genes.
eg.-CO4.Ctl.DGE.csv, here the 2 conditions are CO4(pair1) and Ctl(pair2) and field baseMeanA means average count for CO4(pair1) and baseMeanB means average count for Ctl(pair2).baseMean is average of baseMeanA and baseMeanB.
<pair1.pair2>.DGE.update_annot.csv -Same information as <pair1.pair2>.DGE.csv with an additional Field "Annotation". The Annotation field shows annotation based on Annotations downloaded from PhytozomeV11. Annotations were missing for 7105 genes which were scanned against PfamA dbase to update the annotations.
eg.-CO4.Ctl.DGE.update_annot.csv

<pair1.pair2>.norm.csv-Normalised count values using DESEQ package
eg.-CO4.Ctl.norm.csv

<pair1.pair2>.MA.jpg- log2foldchange is log2(pair2/pair1).plot of log2foldchange vs. mean of normalised counts for each gene is plotted. The genes shown in RED dots are significant(FDR<0.1)
eg.-CO4.Ctl.MA.jpg
<pair1.pair2>.displot.jpg-dispersion(variance see the DESEQ manual) vs. mean of normalised counts is plotted for each gene
eg.-CO4.Ctl.displot.jpg
<pair1.pair2>.pval.jpg-pvalue distribution of genes is shown. usually Most of the genes would have higher pvalue to be termed significant.
eg.-CO4.Ctl.pval.jpg

temp/ directory contains Any Extra/Modified/temporary files for Record keeping. 
