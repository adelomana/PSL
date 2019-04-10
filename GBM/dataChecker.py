import sys,numpy

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

print(len(arrayGenes))

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
            if sampleIDs[i] not in arrayExpression.keys():
                arrayExpression[sampleIDs[i]]={}

            arrayExpression[sampleIDs[i]][geneName]=rowExpression[i]

print(len(rnaseqGenes))

# convert names
intersect=list(set(arrayGenes) & set(rnaseqGenes))
print(len(intersect))

# remove genes not common
            

# check the rank correlation in within and across technologies

# decide to merge and use quantile normalization
