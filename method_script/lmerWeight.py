#Script that generates the weight for each 10mer.
import operator
import gkmer_methods
import sys

'''
Recursive method that generates all possible l-mers of length l
Inputs
lmerList:	List from the previous iteration of the method. Feed in a blank list as
			input on the first iteration.
l:			Number of iterations left until the full list is generated

Output
lmerList:	Full list of l-mers of length l
'''
def generateLmer(lmerList, l):
    NTList = ['T','C','G','A']
    if(k==0):
        return lmerList
    lmerListNew = list()
    for lmer in lmerList:
        for NT in NTList:
            lmerListNew.append(lmer+NT)
    return generateLmer(lmerListNew, k-1)


'''
Uses gapped k-mer weights in order to calculate the weights for all possible l-mers,
given values for k and l.

Inputs
subLength:			Length of the lmers that we are generating weights for.
numGaps:			How many gaps are present for the gapped k-mers in the weights file.
gkmerWeightsDict:	Dictionary that contains the gapped k-mers as keys, and their
					appropriate weights as the values. The keys does not need to
					exhaust all possible gapped k-mers, but only needs to contain the
					ones deemed "important".
verbose:			Boolean that denotes whether the method should print strings.

Output
weightList			List of tuples with (l-mer, weight), where the weight for the l-mer
					is the sum of the gapped k-mer weights that can be created by the
					l-mer. The list is sorted in order of descending weight.
'''
def generateLmerWeights(subLength, numGaps, gkmerWeightsDict, verbose):
    lmerList = generateLmer([""], subLength)
    
    if verbose:
    	print("Calculating l-mer weights.")
    weightList = list()
    i = 1
    for lmer in lmerList:
        if(i%50000 == 0 and verbose):
            print("Currently on lmer: " + str(i))
        weight = 0
        kmerComp = gkmer_methods.getComplement(lmer)
        gkmerSet = set(gkmer_methods.gkmers(lmer, lmerLength, numGaps))
        gkmerSet.update(gkmer_methods.gkmers(lmerComp, lmerLength, numGaps))
        for gkmer in gkmerSet:
            if(gkmer in gkmerWeightsDict):
                weight += gkmerWeightsDict[gkmer]
        weightList.append((lmer,weight))
        i += 1
        
    weightList = sorted(weightList,key = lambda x:x[1],reverse=True)
    return weightList
    