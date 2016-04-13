#encoding=gbk
f=open('sele_features.txt')
DNBs=set([])
while True:
    line=f.readline()
    if len(line)==0:
        break
    fields=line.split('\r')
    t1=fields[0]
    DNBs.add(t1.upper())
f.close()

f=open('mmuAliase.txt')
idx={}
l=0
while True:
    line=f.readline()
    if len(line)==0:
        break
    l+=1
    if l%1000000==0:
        print l
    fields=line.split('\t')
    t1=fields[1]
    t2=fields[2]
    if(t2.upper() in DNBs):
        idx[t1]=t2.upper()
f.close()

f=open('10090.protein.links.v10.txt')
l=0
lw=0
genes=set([])
edges=set([])
f.readline()
while True:
    line=f.readline()
    if len(line)==0:
        break
    l+=1
    if l%1000000==0:
        print l
    fields=line.split(' ')
    t1=fields[0]
    t2=fields[1]
    t1=t1.split('.')
    t2=t2.split('.')
    t1=t1[1]
    t2=t2[1]
    if((t1 in idx) and (t2 in idx)):
        if idx[t1]==idx[t2]:
            continue
        if (idx[t1] in edges) and idx[t2]==edges[idx[t1]]:
            continue
        if (idx[t2] in edges) and idx[t1]==edges[idx[t2]]:
            continue
        lw+=1
        genes.add(idx[t1])
        genes.add(idx[t2])
        edges.add(tuple(sorted((idx[t1],idx[t2]))))
f.close()
print lw,len(genes),len(edges)

fw=open('ppi_network','w')
for e in edges:
    fw.write('%s\t%s\n'%(e[0],e[1]))
fw.close()

fw=open('ppi_genes','w')
for e in genes:
    fw.write('%s\n'%e)
fw.close()
