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
Inputs
subLength:			Length of the lmers that we are generating weights for.
numGaps:			How many gaps are present for the gapped k-mers in the weights file.
weightsFilePath:	String that gives pathing to the file that has the gapped k-mers
					along with their associated weights.
verbose:			Boolean that denotes whether the method should print strings.

Output
weightList			List of tuples with (l-mer, weight), where the weight for the l-mer
					is the sum of the gapped k-mer weights that can be created by the
					l-mer. The list is sorted in order of descending weight.
'''
def generateLmerWeights(subLength, numGaps, weightsFilePath, verbose):
	if verbose:
    	print("Loading in gkmer weights")
    weightsFile = open(weightsFilePath)
    gkmerWeightsDict = {}
    for line in weightsFile:
        splitLine = line.strip().split()
        gkmer = splitLine[0]
        weight = eval(splitLine[1])
        gkmerWeightsDict[gkmer] = weight
    weightsFile.close()
    if verbose:
    	print('Number of gkmers loaded in: ' + str(len(gkmerWeightsDict)))
    
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
    