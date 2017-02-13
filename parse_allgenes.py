#!/usr/bin/env python
import re
import sys
fh=open(sys.argv[1])
fh2=open(sys.argv[2],'w')
header=fh.readline()
print header
fh2.write(header)
for line in fh:
    line=line.strip('\n').strip()
    gid_def=line.split('\t')
    #print len(gid_def)
    gid=gid_def[0].strip()
    print gid
    des=gid_def[1:]
    print des
    #secondField=des[0].strip()
    if des==[]:
        print des
        print gid
        fh2.write(gid+'\t'+'NA'+'\t'+'NA'+'\n')
    else:
        fh2.write(gid+'\t'+des[0].strip('\t').strip()+'\t'+des[1].strip('\t').strip()+'\n')
fh2.close()
fh.close()
