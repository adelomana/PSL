import sys,os
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

import miner2
import miner2.preprocess

# 0. user defined variables
expressionFile='/Volumes/omics4tb2/alomana/projects/PSL/MM/data/IA12Zscore.csv'
resultsDir='/Volumes/omics4tb2/alomana/projects/PSL/MM/results/'

# STEP 0: load the data
expressionData, conversionTable = miner2.preprocess.main(expressionFile)

individual_expression_data = [expressionData.iloc[:,i] for i in range(50)]
matplotlib.pyplot.boxplot(individual_expression_data)
matplotlib.pyplot.title("Patient expression profiles",FontSize=14)
matplotlib.pyplot.ylabel("Relative expression",FontSize=14)
matplotlib.pyplot.xlabel("Sample ID",FontSize=14)

figureName=resultsDir+'figure.xx.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()

print('figrue done')s
#plt.figure()
#_ = plt.hist(expressionData.iloc[0,:],bins=100,alpha=0.75)
#plt.title("Expression of single gene",FontSize=14)
#plt.ylabel("Frequency",FontSize=14)
#plt.xlabel("Relative expression",FontSize=14)
#plt.figure()
#_ = plt.hist(expressionData.iloc[:,0],bins=200,color=[0,0.4,0.8],alpha=0.75)
#plt.ylim(0,350)
#plt.title("Expression of single patient sample",FontSize=14)
#plt.ylabel("Frequency",FontSize=14)
#plt.xlabel("Relative expression",FontSize=14)

"""

# 1. step 1: clustering
#step1clustering()

