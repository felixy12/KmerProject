import numpy as np
import gkmerHelpers
import time

'''
Uses the dictionary of important words, along with a new fasta file in order to 
create a matrix, where the rows denote sequences, and columns denote the feature counts.

Inputs
subLength:		Length of the gapped k-mer, l
numGap:			Number of gaps in the gapped k-mer. k = subLength - numGap
seqFilePath:	Path to the fasta file that contains the sequences we are generating
				feature vectors for.
importantWords:	Dictionary that contains gapped k-mers that exceed a certain zScore, and
				their index within the feature vector.
verbose:		Boolean flag that determins whether certain strings are printed.
'''
def getFeatures(subLength, numGap, seqFilePath, importantWords, verbose):
    weightsFile = open(weightsFilePath, 'r')

    seqNameList, seqList = gkmer_methods.loadFasta(seqFilePath)
    importantWords = {}

    start = time.time()
    print("Creating feature matrix.")
    ftMat = np.zeros((len(seqList), len(importantWords)), dtype = 'int8')
    for seqInd in range(len(seqList)):
    	curSeq = seqInd+1
        ftMat[seqInd,:] = gkmer_methods.featureVector(seqList[seqInd], subLength, numGap, importantWords)
        if(curSeq%500 == 0 and verbose):
            print("Currently processing sequence: " + str(seqInd))
    end = time.time()
    print("Done. Time Elapsed: " + str(int((end-start)/60)) + " minutes " + str(int((end-start)%60)) + " seconds") 
    return ftMat