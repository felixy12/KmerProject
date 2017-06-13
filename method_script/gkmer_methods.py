import itertools
import sys
import numpy as np
import math

'''
Loads in data from a .fa file.

Inputs
fileName:	String containing the name of the fasta file

Outputs
seqNameList:	A python List of the names of the sequences in the fasta file, in order
				of how they are in the file.
seqList:		A python List of the sequences themselves in the fasta file, in order
				of how they are in the file. This means that a sequence and it's name
				are related by what index they are in the lists.
'''
def loadFasta(fileName):
    seqFile = open(fileName, 'r')
    seqNameList = list()
    seqList = list()
    seq = seqFile.readline()
    while(">" not in seq and not seq==""):
        seq = seqFile.readline()
    seqNameList.append(seq.strip()[1:])
    seqLine = ""
    for seq in seqFile:
        if(">" not in seq and not seq==""): 
            seqLine = seqLine+seq.strip()
        elif(">" in seq):
            seqNameList.append(seq.strip()[1:])
            seqLine = seqLine.upper()
            seqList.append(seqLine)
            seqLine = ""
    seqLine = seqLine.upper()
    seqList.append(seqLine)
    seqFile.close()
    return seqNameList, seqList

'''
Loads in support vector alpha values from file. 

Inputs
fileName:	String that contains the name of the svalpha file

Outputs
alphaList: 	A python List that contains all alpha values in order in which they appeared
			in the file. When run with loadFasta on svseq.fa, the indices in alphaList
			will correspond to the indices outputted by loadFasta
'''
#Given a string that points to the svalpha file, it'll load in all alpha
#values in a list in order of the file.
def loadAlpha(fileName):
    alphaFile = open(fileName, 'r')
    alphaList = list()
    for aLine in alphaFile:
        aLine = aLine.strip().split()
        alpha = float(aLine[1])
        alphaList.append(alpha)
    alphaFile.close()
    return alphaList

'''
Creates the reverse complement of a given DNA sequence. If a gapped sequence is given,
the gaps are ignored (they stay as gaps).

Inputs
seq:		A string containing the DNA sequence

Outputs
seq:		A string containing the reverse complement
'''
def complement(seq):
    seq = seq[::-1]
    seq = seq.replace('A', 'a')
    seq = seq.replace('T', 'A')
    seq = seq.replace('a', 'T')
    seq = seq.replace('C', 'c')
    seq = seq.replace('G', 'C')
    seq = seq.replace('c', 'G')
    return seq

'''
Generates a list of gapped k-mers from a given sequence.

Inputs
seq:		String that contains the sequence to be processed into gkmers
subLength:	Length of the l-mers that the sequence will be broken down into
numGap:		Number of gaps that are present within each l-mer. A gapped k-mer
			then has k = subLength-numGap
			
Outputs
gappedList:	A list containing all of the gapped k-mers that can be created
			by the sequence. This list includes repeats.
'''
def gkmers(seq, subLength, numGap):
    subList = list()
    for i in range(0, len(seq)-subLength+1):
    	subList.append(seq[i: i+subLength])
    comb = itertools.combinations(range(0,subLength), numGap)
    indList = list(comb)
    gappedList = list()
    for subS in subList:
        for iComb in indList:
        	tempSubS = subS
        	for i in iComb:
        		tempSubS = tempSubS[:i]+'-'+tempSubS[i+1:]
            gappedList.append(tempSubS)
    return gappedList

'''
Calculates the contribution from a given sequence to the weight of the gapped k-mers.

Inputs
seq:		String that contains the sequence to be processed into gkmers
subLength:	Length of the l-mers that the sequence will be broken down into
numGap:		Number of gaps that are present within each l-mer. A gapped k-mer
			then has k = subLength-numGap
alpha:		The support vector alpha value of the sequence.
		
Outputs
countDict:	Python dictionary with the gapped k-mer as the key, and the contribution
			of the inputted sequence on the weight of the gapped k-mer as the value.
			The contribution for any one gk-mer is (num_count/norm_factor)*alpha
'''
def gkmerWeightFromSeq(seq, subLength, numGap, alpha):
    subList = list()
    for i in range(0, len(seq)-subLength+1):
    	subList.append(seq[i: i+subLength])
    comb = itertools.combinations(range(0,subLength), numGap)
    indList = list(comb)
    countDict = {}
    for subS in subList:
        for iComb in indList:
        	gkmer = subS
        	for i in iComb:
        		gkmer = gkmer[:i]+'-'+gkmer[i+1:]
            countDict[gkmer]=countDict.get(kmer,0.0)+alpha
            
    Values = np.square(countDict.values())    
    normFactor = np.sqrt(np.sum(Values))
    for key in countDict.keys():
        countDict[key] = (countDict[key]/normFactor)*alpha
    return countDict

'''
Counts the number of times an important gapped k-mer appears in a sequence. Essentially
creates the feature vector for a sequence, where each feature is the number of times
a certain gapped k-mer appears in the sequence.

Inputs
seq:			String that contains the sequence to be processed into gkmers
subLength:		Length of the l-mers that the sequence will be broken down into
numGap:			Number of gaps that are present within each l-mer. A gapped k-mer
				then has k = subLength-numGap
importantWords:	A dictionary of gapped k-mers that are deemed "important" 
				(has z-score above a certain threshold). The keys to the dictionary 
				are the important gapped k-mers, and the values are the indices
				associated with them.
'''
def getImportantCounts(seq, subLength, numGap, importantWords):
	counts = np.zeros(len(importantWords), dtype = 'int8')
	subList = list()
    for i in range(0, len(seq)-subLength+1):
    	subList.append(seq[i: i+subLength])
    comb = itertools.combinations(range(0,subLength), numGap)
    indList = list(comb)
    countDict = {}
    for subS in subList:
        for iComb in indList:
        	gkmer = subS
        	for i in iComb:
        		gkmer = gkmer[:i]+'-'+gkmer[i+1:]
        		if(gkmer in importantWords):
        			index = importantWords.get(gkmer)
        			counts[index] = counts[index]+1
    return counts