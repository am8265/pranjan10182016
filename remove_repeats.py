#!/usr/bin/env python
import re
import sys
fh=open(sys.argv[1])
fh2=open(sys.argv[2],'w')
header=fh.readline()
fh2.write(header)
geneDescription=dict()
for line in fh:
    line=line.strip('\n').strip()
    gene_des=line.split('\t')
    gene=gene_des[0].strip()
    des=gene_des[1].strip()
    if geneDescription.get(gene,'0')=='0':#i.e the geneID Never Existed!
        geneDescription[gene]=[des]
    elif geneDescription.get(gene)==['NA'] and des!='NA':#i.e presently we have a Valid Annotation!
        geneDescription[gene]=[des]
    elif geneDescription.get(gene)!=['NA'] and des=='NA':#i.e previous annotation was better 
        pass
    elif geneDescription.get(gene)==[des]:#2 situations where both 1) are NA 2) or same repeatin
        pass
    else:
        geneDescription[gene].append(des)
        geneDescription[gene]=list(set(geneDescription[gene]))
fh.close()
for gene in geneDescription.keys():
    des=''
    if geneDescription[gene]>1:#i.e that id has more than on e definition
        for descrip in geneDescription[gene]:
            des=des+descrip
    else:
        des=geneDescription[gene][0]
    fh2.write(gene+'\t'+des+'\t'+des+'\n')
fh2.close()
