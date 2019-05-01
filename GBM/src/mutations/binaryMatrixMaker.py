###
### This script reads TCGA mutation files from TCGA and creates a binary matrix of patients and genes for mutation presence.
###

import os,sys

def TCGAfileReader(path,mutationTypes):

    '''
    This function reads TCGA mutation files and provides a list of mutations
    '''

    mutationSet=[]

    patientLabel=path.split('/')[-1].split('.maf')[0]
    #print('\t working with {}'.format(patientLabel))
    
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
TCGAfolder='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/gdac.broadinstitute.org_GBM.Mutation_Packager_Calls.Level_3.2016012800.0.0/'
resultsFile='/Volumes/omics4tb2/alomana/projects/PSL/GBM/data/mutations/mutations.GBM.TCGA.2019.04.30.csv'

# 1. read data and create list of tuples with mutations
# mutations=[(patient,gene), (), ...]
print('read data...')

mutations=[]
mutationTypes=[]
patientLabels=[]

detectedFiles=os.listdir(TCGAfolder)
patientFiles=[element for element in detectedFiles if 'TCGA' in element and '._TCGA' not in element]

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
mutatedGeneList.remove('Unknown')
for geneLabel in mutatedGeneList[:5]:
    print(geneLabel,mutatedGenes[geneLabel])

# 2. create binary csv file
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
