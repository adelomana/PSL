

with open('/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/mutations/mutations.GBM.TCGA.2019.05.01.csv','r') as f:
    line=f.readline()
    v=line.split(',')
    v[-1]=v[-1].replace('\n','')
    v.pop(0)
    print(len(v))

with open('/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/normalization/QN.MA.data.csv','r') as f:
    line=f.readline()
    x=line.split(',')
    x[-1]=x[-1].replace('\n','')
    x.pop(0)
    print(len(x))


i=list(set(x) & set(v))
print(len(i))
