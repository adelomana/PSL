import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

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
    
    # incorporate new values based on position
    normalized={}
    for sample in samples:
        normalized[sample]={}
        values=[]
        for gene in genes:
           values.append(adata[sample][gene])        
        order=numpy.argsort(values)
        ranks = order.argsort()
        newValues=averages[ranks]
        for i in range(len(genes)):
            normalized[sample][genes[i]]=newValues[i]
        
    return normalized

def histogrammer(theData):

    '''
    This function creates a histogram.
    '''    

    x=[]; y=[]
    
    n,bins=numpy.histogram(theData,bins=int(numpy.sqrt(len(theData))))

    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        x.append(center)
    x.pop()
    y=numpy.array(n)
    y=list(y/float(sum(y)))

    return x,y

# 0. user defined variables
# microarray data file
arrayDataFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/microarray/GBM.medianexp.txt'

# RNA-seq, normalized counts
rnaseqDataFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/RNAseq/GBM.uncv2.mRNAseq_RSEM_all.txt'

# 1. read data
print('reading data...')

arrayExpression={}; rnaseqExpression={}

# 1.1. read array data
arrayGenes=[]
with open(arrayDataFile,'r') as f:

    headerLine=f.readline()
    v=headerLine.split('\t')
    sampleIDs=v[1:]
    sampleIDs[-1]=sampleIDs[-1].replace('\n','')
    
    next(f)

    for line in f:
        v=line.split('\t')
        geneName=v[0]

        if geneName not in arrayGenes:
            arrayGenes.append(geneName)
            
        rowExpression=[float(element) for element in v[1:]]

        # create dictionary
        for i in range(len(sampleIDs)):
            if sampleIDs[i] not in arrayExpression.keys():
                arrayExpression[sampleIDs[i]]={}

            arrayExpression[sampleIDs[i]][geneName]=rowExpression[i]

# 1.2. read RNAseq data
rnaseqGenes=[]
with open(rnaseqDataFile,'r') as f:

    headerLine=f.readline()
    v=headerLine.split('\t')
    sampleIDs=v[1:]
    sampleIDs[-1]=sampleIDs[-1].replace('\n','')

    for line in f:
        v=line.split('\t')
        geneName=v[0].split('|')[0]
        if geneName not in rnaseqGenes:
            rnaseqGenes.append(geneName)
        rowExpression=[float(element) for element in v[1:]]
        log2plusOne=[numpy.log2(element+1) for element in rowExpression]
        
        # create dictionary
        for i in range(len(sampleIDs)):
            if sampleIDs[i] not in rnaseqExpression.keys():
                rnaseqExpression[sampleIDs[i]]={}

            rnaseqExpression[sampleIDs[i]][geneName]=log2plusOne[i]

# convert names
commonGenes=list(set(arrayGenes) & set(rnaseqGenes))
commonGenes.sort()

# remove genes not common in array
genes2remove=[]
oneSample=list(arrayExpression.keys())[0]

for gene in arrayExpression[oneSample]:
    if gene not in commonGenes:
        genes2remove.append(gene)
for sample in arrayExpression:
    for gene in genes2remove:
        del arrayExpression[sample][gene]

# remove genes not common in array
genes2remove=[]
oneSample=list(rnaseqExpression.keys())[0]

for gene in rnaseqExpression[oneSample]:
    if gene not in commonGenes:
        genes2remove.append(gene)
for sample in rnaseqExpression:
    for gene in genes2remove:
        del rnaseqExpression[sample][gene]


# 2. analysis
print()
print('analysis...')

fewSamplesArray=list(arrayExpression.keys())[:2000]
fewSamplesRNAseq=list(rnaseqExpression.keys())[:2000]
fewGenes=commonGenes[:2000]

"""

### ARRAY
print('\t Arrays...')

# get the mean and std across genes
values=[]
for gene in fewGenes:
    for sample in fewSamplesArray:
        values.append(arrayExpression[sample][gene])
    low=numpy.min(values)
    high=numpy.max(values)
    mean=numpy.mean(values)
    stdDev=numpy.std(values)
    variance=numpy.var(values)

# get the mean and std across samples
values=[]
for sample in fewSamplesArray:
    for gene in fewGenes:
        values.append(arrayExpression[sample][gene])
    low=numpy.min(values)
    high=numpy.max(values)
    mean=numpy.mean(values)
    stdDev=numpy.std(values)
    variance=numpy.var(values)

### RNAseq
print('\t RNAseq...')

# get the mean and std across genes
values=[]
for gene in fewGenes:
    for sample in fewSamplesRNAseq:
        values.append(rnaseqExpression[sample][gene])
    low=numpy.min(values)
    high=numpy.max(values)
    mean=numpy.mean(values)
    stdDev=numpy.std(values)
    variance=numpy.var(values)

# get the mean and std across samples
values=[]
for sample in fewSamplesRNAseq:
    for gene in fewGenes:
        values.append(rnaseqExpression[sample][gene])
    low=numpy.min(values)
    high=numpy.max(values)
    mean=numpy.mean(values)
    stdDev=numpy.std(values)
    variance=numpy.var(values)

"""

# 2.1 merge two data sets
fullSet={}
for sample in fewSamplesArray:
    fullSet[sample]={}
    for gene in fewGenes:
        fullSet[sample][gene]=arrayExpression[sample][gene]
for sample in fewSamplesRNAseq:
    fullSet[sample]={}
    for gene in fewGenes:
        fullSet[sample][gene]=rnaseqExpression[sample][gene]
        
# 2.2. plot distributions
for sample in fewSamplesArray:
    values=[]
    for gene in fewGenes:
        values.append(arrayExpression[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='blue',lw=1,alpha=1/2)

for sample in fewSamplesRNAseq:
    values=[]
    for gene in fewGenes:
        values.append(rnaseqExpression[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='red',lw=1,alpha=1/2)

# 2.2. run quantile normalization
print('quantile normalization...')
quantileNormalized=expressionQuantileNormalization(fullSet)
print('... completed.')

print('plotting...')
for sample in quantileNormalized:
    values=[]
    for gene in quantileNormalized[sample]:
        values.append(quantileNormalized[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='black',lw=1,alpha=1/2)
    
matplotlib.pyplot.xlabel('expression')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.expression.distribution.pdf')
matplotlib.pyplot.clf()
    
# 3. store data into csv
print('storing...')

# 3.1. store full data set
sampleNames=list(quantileNormalized.keys())
sampleNames.sort()
geneNames=list(quantileNormalized[sampleNames[0]].keys())
geneNames.sort()

fileName='results/full.QN.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(quantileNormalized[sampleName][geneName]))
        
        f.write('\n')

# 3.2. store MA
sampleNames=list(arrayExpression.keys())
sampleNames.sort()
geneNames=list(arrayExpression[sampleNames[0]].keys())
geneNames.sort()

fileName='results/MA.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(arrayExpression[sampleName][geneName]))
        
        f.write('\n')

# 3.3. store RS
sampleNames=list(rnaseqExpression.keys())
sampleNames.sort()
geneNames=list(rnaseqExpression[sampleNames[0]].keys())
geneNames.sort()

fileName='results/RS.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(rnaseqExpression[sampleName][geneName]))
        
        f.write('\n')
        

