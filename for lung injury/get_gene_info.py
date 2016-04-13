import random
f=open('new_ppi_network2')
gene_info={}
edges=[]
while True:
    line=f.readline()
    if len(line)==0:
        break
    fields=line.split('\t')
    g1=fields[0]
    g2=fields[1]
    edges.append((g1,g2))
    if g1 not in gene_info:
        gene_info[g1]=0
    tmp=gene_info[g1]
    tmp=tmp+1
    gene_info[g1]=tmp
    if g2 not in gene_info:
        gene_info[g2]=0
    tmp=gene_info[g2]
    tmp=tmp+1
    gene_info[g2]=tmp
f.close()

gene_info=sorted(gene_info.iteritems(),key=lambda d:d[1], reverse = True)

fw=open('gene_flag.txt','w')
DNB=set([])
for item in gene_info:
    gene=item[0]
    count=item[1]
    if count>=6:
        flag=1
        state1=random.uniform(1, 10)
        state2=random.uniform(1, 10)
        state3=random.uniform(30, 40)
        if count>=11:
            flag=2
            state1=random.uniform(30, 40)
            state2=random.uniform(30, 40)
            state3=random.uniform(1, 10)
        DNB.add(gene)
    else:
        flag=0
        state1=random.uniform(1, 20)
        state2=random.uniform(1, 20)
        state3=random.uniform(1, 20)
    fw.write('%s\t%f\t%f\t%f\t%d\t%d\n'%(gene,state1,state2,state3,count,flag))
fw.close()

#fw=open('network_info.txt','w')
for e in edges:
    g1=e[0]
    g2=e[1]
    if g1 in DNB and g2 in DNB:
        pcc1=random.uniform(0.1, 0.3)
        pcc2=random.uniform(0.9, 1)
        pcc3=random.uniform(0.1, 0.3)
    elif g1 not in DNB and g2 not in DNB:
        pcc1=random.uniform(0.1, 0.3)
        pcc2=random.uniform(0.1, 0.3)
        pcc3=random.uniform(0.1, 0.3)
    else:
        pcc1=random.uniform(0.3, 0.6)
        pcc2=random.uniform(0.1, 0.3)
        pcc3=random.uniform(0.3, 0.6)
        
    #fw.write('%s (pp) %s\t%f\t%f\t%f\n'%(g1,g2,pcc1,pcc2,pcc3))
#fw.close()

