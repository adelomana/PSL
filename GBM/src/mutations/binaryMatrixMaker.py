###
### This script reads TCGA mutation files from TCGA and creates a binary matrix of patients and genes for mutation presence.
###

import os,sys

def TCGAfileReader(path,mutationTypes):

    '''
    This function reads TCGA mutation files and provides a list of mutations
    '''

    mutationSet=[]

    patientLabel='-'.join(path.split('/')[-1].split('-')[:3])
    
    with open(path,'r') as f:
        next(f)
        for line in f:
            v=line.split('\t')

            if len(v) > 10:

                geneLabel=v[0]
                mutationType=v[8]
            
                if mutationType not in mutationTypes:
                    mutationTypes.append(mutationType)

                if mutationType != 'Silent':
                    mutationSet.append((patientLabel,geneLabel))

    return mutationSet,mutationTypes,patientLabel

# 0. user defined variables
TCGAfolder='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/mutations/gdac.broadinstitute.org_GBM.Mutation_Packager_Calls.Level_3.2016012800.0.0/'
resultsFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/mutations/mutations.GBM.TCGA.2019.05.01.csv'
expressionDataFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/normalization/QN.MA.data.csv'

# 1. read data and create list of tuples with mutations
# mutations=[(patient,gene), (), ...]
print('read data...')

mutations=[]
mutationTypes=[]
patientLabels=[]

# 1.1. determine mutation patient files
detectedFiles=os.listdir(TCGAfolder)
allPatientFiles=[element for element in detectedFiles if 'TCGA' in element and '._TCGA' not in element]

# 1.2. drop patients for which no expression data is available
with open(expressionDataFile,'r') as f:
    line=f.readline()
    v=line.split(',')
    expressionPatientIDs=v[1:]
    expressionPatientIDs[-1]=expressionPatientIDs[-1].replace('\n','')

patientFiles=[]
for putative in allPatientFiles:
    ID='-'.join(putative.split('-')[:3])
    if ID in expressionPatientIDs:
        patientFiles.append(putative)

# 1.3. read info in mutation patient files
for patientFile in patientFiles:
    path=TCGAfolder+patientFile
    mutationSet,mutationTypes,patientLabel=TCGAfileReader(path,mutationTypes)
    for element in mutationSet:
        if element not in mutations:
            mutations.append(element)
    patientLabels.append(patientLabel)
    
mutationTypes.sort()
print(mutationTypes)

mutatedGenes={}
for item in mutations:
        if item[1] in mutatedGenes:
                mutatedGenes[item[1]]=mutatedGenes[item[1]]+1
        if item[1] not in mutatedGenes:
                mutatedGenes[item[1]]=1

mutatedGeneList=[k for k in sorted(mutatedGenes, key=mutatedGenes.get, reverse=True)]
mutatedGeneList=[element for element in mutatedGeneList if mutatedGenes[element] >= 3] # the mutation needs to be present in at least 3 patients to be considered
mutatedGeneList.remove('Unknown')

for i in range(len(mutatedGeneList)):
    print(i,mutatedGeneList[i],mutatedGenes[mutatedGeneList[i]])

# 3. create binary csv file
print('Store binary csv file...')
with open(resultsFile,'w') as f:

    for patientLabel in patientLabels:
        f.write(',{}'.format(patientLabel))
    f.write('\n')

    for geneLabel in mutatedGeneList:
        f.write(geneLabel)

        for patientLabel in patientLabels:
            if (patientLabel,geneLabel) in mutations:
                f.write(',1')
            else:
                f.write(',0')
        f.write('\n')
