###
### This script runs MINER for GBM data from TCGA
###

import sys

# 0. preliminaries

# 0.1. add MINER source to Python path
sys.path.append('/Users/alomana/github/miner.MW/miner/src')
import miner

# 0.2. import packages required by MINER
import os,pandas,json,seaborn,time,pickle,collections
import scipy,scipy.stats
import numpy,numpy.random
import sklearn,sklearn.decomposition,sklearn.manifold
import multiprocessing,multiprocessing.pool
import matplotlib,matplotlib.pyplot

# 0.3. user defined variables
inputDataFile='/Users/alomana/github/PSL/GBM/results/MA.data.csv'
outputDir='/Volumes/omics4tb2/alomana/projects/PSL/GBM/results/MINER/'
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

# 1. input and preprocess expression data
expressionData=pandas.read_csv(inputDataFile,header=0,index_col=0,sep=",")
expressionData,conversionTable=miner.identifierConversion(expressionData)

