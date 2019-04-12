import numpy,sys


def expressionQuantileNormalization(adata):


    '''
    Quantile normalization as in wikipedia
    '''
    
    samples=list(adata.keys())
    genes=list(adata[samples[0]].keys())
    
    # rank genes within each sample and compute their mean
    sortedMatrix=[]
    for sample in samples:
        values=[]
        for gene in genes:
            value=adata[sample][gene]
            values.append(value)
        values.sort()
        sortedMatrix.append(values)
    averages=numpy.mean(numpy.array(sortedMatrix),axis=0)
    print(averages)
    
    # incorporate new values based on position
    normalized={}
    for sample in samples:
        normalized[sample]={}
        values=[]
        for gene in genes:
           values.append(adata[sample][gene])        
        order=numpy.argsort(values)
        ranks = order.argsort()
        print(sample,values,order,ranks)
        newValues=averages[ranks]
        print('new',newValues)
        for i in range(len(genes)):
            normalized[sample][genes[i]]=newValues[i]
        
    return normalized

adata={}
adata['sample1']={}
adata['sample2']={}
adata['sample3']={}
adata['sample4']={}

adata['sample1']['geneA']=2
adata['sample1']['geneB']=5
adata['sample1']['geneC']=4
adata['sample1']['geneD']=3
adata['sample1']['geneE']=3

adata['sample2']['geneA']=4
adata['sample2']['geneB']=14
adata['sample2']['geneC']=8
adata['sample2']['geneD']=8
adata['sample2']['geneE']=9

adata['sample3']['geneA']=4
adata['sample3']['geneB']=4
adata['sample3']['geneC']=6
adata['sample3']['geneD']=5
adata['sample3']['geneE']=3

adata['sample4']['geneA']=5
adata['sample4']['geneB']=7
adata['sample4']['geneC']=9
adata['sample4']['geneD']=8
adata['sample4']['geneE']=5

print(adata)
quantileNormalized=expressionQuantileNormalization(adata)
print(quantileNormalized)
