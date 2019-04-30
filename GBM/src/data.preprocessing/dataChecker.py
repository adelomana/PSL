import sys,numpy
import sklearn,sklearn.preprocessing
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

def expressionScaling(adata):

    '''
    Bring data to a range betwenn 0 and 1.
    '''

    samples=list(adata.keys())
    genes=list(adata[samples[0]].keys())
    
    # find min and max
    allValues=[]
    for sample in samples:
        for gene in genes:
            value=adata[sample][gene]
            allValues.append(value)
    low=min(allValues)
    high=max(allValues)

    # create a new dictionary
    scaled={}
    for sample in samples:
        scaled[sample]={}
        for gene in genes:
            value=adata[sample][gene]
            newValue=(value-low)/(high-low)
            scaled[sample][gene]=newValue

    # find min and max
    allValues=[]
    for sample in samples:
        for gene in genes:
            value=scaled[sample][gene]
            allValues.append(value)
    low=min(allValues)
    high=max(allValues)
    print('after',low,high)
            
    return scaled
                
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

###
### MAIN
###

### 0. user defined variable
### 1. read data
### 2. normalizations
### 3. plot original and nomalized distributions
### 4. visualize samples umap
### 4. visualize samples heatmap
### 4. store data

# 0. user defined variables
# microarray data file
arrayDataFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/microarray/GBM.medianexp.txt'

# RNA-seq, normalized counts
rnaseqDataFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/RNAseq/GBM.uncv2.mRNAseq_RSEM_all.txt'

# results folder
resultsFolder='/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/normalization/'

# 1. read data
print('Reading data...')

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

# 1.3. find common intersect
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

# 2. normalizations
print('Quantile normalization...')

MAsamples=list(arrayExpression.keys())[:]
RSsamples=list(rnaseqExpression.keys())[:]
geneNames=commonGenes[:]

# 2.1 merge two data sets
fullSet={}
for sample in MAsamples:
    fullSet[sample]={}
    for gene in geneNames:
        fullSet[sample][gene]=arrayExpression[sample][gene]
for sample in RSsamples:
    fullSet[sample]={}
    for gene in geneNames:
        fullSet[sample][gene]=rnaseqExpression[sample][gene]
        
# 2.2. transform array data into matrix
allValues=[]
for sample in MAsamples:
    values=[]
    for gene in geneNames:
        value=arrayExpression[sample][gene]
        values.append(value)
    allValues.append(values)
x=numpy.array(allValues)
arrayExpressionNum=numpy.transpose(x)

# 2.3. run quantile normalization
QNF=expressionQuantileNormalization(fullSet)
QNMA=expressionQuantileNormalization(arrayExpression)
QNRS=expressionQuantileNormalization(rnaseqExpression)

QNMAscaled=expressionScaling(QNMA)

# 2.4. run sklearn quantile normalization
sklearn.preprocessing.quantile_transform(arrayExpressionNum,axis=0,output_distribution='normal')

print('\t ... completed.')

# 3. plot original and normalized distributions
print('Plotting...')

# 3.1. plot original distributions
for sample in MAsamples:
    values=[]
    for gene in geneNames:
        values.append(arrayExpression[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='blue',lw=1,alpha=0.1)

for sample in RSsamples:
    values=[]
    for gene in geneNames:
        values.append(rnaseqExpression[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='red',lw=1,alpha=0.1)

for sample in QNMA:
    values=[]
    for gene in QNMA[sample]:
        values.append(QNMA[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='cyan',lw=2,alpha=0.1)

for sample in QNRS:
    values=[]
    for gene in QNRS[sample]:
        values.append(QNRS[sample][gene])
        
    x,y=histogrammer(values)
    matplotlib.pyplot.plot(x,y,'-',color='magenta',lw=2,alpha=0.1)

matplotlib.pyplot.xlabel('Expression')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(resultsFolder+'figure.expression.original.distributions.pdf')
matplotlib.pyplot.clf()

# 3.2. plot sklearn distributions
for sample in arrayExpressionNum.T:
    x,y=histogrammer(sample)
    matplotlib.pyplot.plot(x,y,'-',color='cyan',lw=1,alpha=0.1)

matplotlib.pyplot.xlabel('Expression')
matplotlib.pyplot.ylabel('Probability')
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(resultsFolder+'figure.expression.sklearn.pdf')
matplotlib.pyplot.clf()

# 4. visualize samples



sys.exit()


# 4. store data into csv
print('storing...')

# 3.1. store original data sets
# store MA
sampleNames=list(arrayExpression.keys())
sampleNames.sort()
geneNames=list(arrayExpression[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'MA.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(arrayExpression[sampleName][geneName]))
        
        f.write('\n')

# store RS
sampleNames=list(rnaseqExpression.keys())
sampleNames.sort()
geneNames=list(rnaseqExpression[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'RS.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(rnaseqExpression[sampleName][geneName]))
        
        f.write('\n')

# store full
sampleNames=list(fullSet.keys())
sampleNames.sort()
geneNames=list(fullSet[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'full.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(fullSet[sampleName][geneName]))
        
        f.write('\n')

# 3.2. store QN data sets
# store MA
sampleNames=list(QNMA.keys())
sampleNames.sort()
geneNames=list(QNMA[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'QN.MA.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(QNMA[sampleName][geneName]))
        
        f.write('\n')

# store RS
sampleNames=list(QNRS.keys())
sampleNames.sort()
geneNames=list(QNRS[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'QN.RS.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(QNRS[sampleName][geneName]))
        
        f.write('\n')

# store full
sampleNames=list(QNF.keys())
sampleNames.sort()
geneNames=list(QNF[sampleNames[0]].keys())
geneNames.sort()

fileName=resultsFolder+'QN.full.data.csv'
with open(fileName,'w') as f:
    header='geneID,'+','.join(sampleNames)
    f.write(header)
    f.write('\n')
    
    for geneName in geneNames:
        f.write(geneName)

        for sampleName in sampleNames:
            f.write(',{}'.format(QNF[sampleName][geneName]))
        
        f.write('\n')
