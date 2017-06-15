import numpy
import sys
import gkmer_methods

'''
Given a fasta file, counts how many times each l-mer appears within the fasta file.

Inputs
subLength:		Length of the subsequence, l.
fastaFilePath:	Path to the fasta file to be counted.
lmerFilePath:	Path to the list of lmers, ranked in order from highest to lowest weight.

Outputs
lmerCountList:	List of tuples, (l-mer, count). The list is in order by highest l-mer
				weight to lowest.
'''
def lmerCounts(subLength, fastaFilePath, lmerFilePath):

    print('Loading in sequences.')
    seqNameList, seqList = gkmerHelpers.loadFasta(fastaFilePath)
    print('Loading in k-mer weights.')
    lmerFile = open(lmerFilePath)
    
    lmerCountDict = {}
    lmerOrder = list()
    for line in lmerFile:
        splitLine = line.strip().split()
        lmerOrder.append(splitLine[0])
        lmerCountDict[splitLine[0]] = 0
    for i in range(len(seqList)):
        if((i+1) % 1000 == 0): 
            print('Currently processing sequence: ' + str(i+1))
        seq = seqList[i]
        lmerList = list()
        for i in range(len(seq)-subLength+1):
        	lmerList.append(seq[i:i+subLength])
        for l in lmerList:
            lmerCountDict[l]+=1
    lmerCountList = list()
    for lmer in lmerOrder:
		lmerCountList.append(lmer, lmerCountDict[lmer])
	return lmerCountList
