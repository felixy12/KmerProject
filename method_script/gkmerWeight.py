import gkmer_methods
import time
import numpy as np

'''
Calculates the weights for all possible gapped k-mers that appear in the SVSeq fasta
file. Reverse complements are combined, where the k-mer that appears 
first alphabetically ('-' occurs before 'A') is the one that is kept.
Inputs
subLength:		Also known as l, denotes how long each gapped k-mer will be
numGap:			Denotes how many gaps there are per k-mer. k = subLength - numGap
seqFilePath:	Path to the directory that contains the sequences for the support
				vectors from the SVM.
alphaFilePath:	Path to the directory that contains the alpha files for the sequences
				in the support vectors of SVM.	
verbose:		Flag on whether certain print line statements will be used

outputs
weightDict		Python dictionary with the gapped k-mer as the key, and the associated
				weight calculated by the support vector sequences as the value.		
'''
def gkmerWeight(subLength, numGap, seqFilePath, alphaFilePath, verbose):
        
    strLength = 0;
    dictSize = 0;
    totalStart = time.time()
    start = totalStart
    end = start
    weightDict = {}
    numPerAnalysis = 500 

    print("Calculating gapped k-mer weights")
    print("Loading files.")
    seqNameList, seqList = gkmer_methods.loadFasta(seqFilePath)
    alphaList = gkmer_methods.loadAlpha(alphaFilePath)
    print("Found " + str(len(seqNameList)) + " sequences to process.")
    print("Generating word weights.")
    for i in range(len(alphaList)):
        currSeq = i+1
        alpha = alphaList[i]
        seqLine = seqList[i]
        tempDict = gkmer_methods.gkmerWeightFromSeq(seqLine, subLength, numGap, alpha)
        dictSize = dictSize + len(tempDict)
        for key in tempDict.keys():
            weightDict[key] = weightDict.get(key, 0) + tempDict[key]
        strLength = strLength+len(seqLine)
        if(currSeq%numPerAnalysis == 0): 
            end = time.time()
            if verbose:
            	print("Summary of sequences " + str(currSeq-numPerAnalysis+1) + "-"+str(currSeq))    
            	print("Time Elapsed: " + str(end-start) + " seconds")
            	print("Average length of sequence: " + str(strLength/numPerAnalysis)+"bp")
            	print("Average size of dictionary for each sequence: " + str(dictSize/numPerAnalysis) + " entries")
            	print("Current number of unique words: " + str(len(weightDict.keys())))
            	print("")
            strLength = 0
            dictSize = 0;
            start = end
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(alphaList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds")

	#Removes any gapped k-mers with an invalid nucleotide, and combines reverse complement weights.
    Keys = weightDict.keys()
    Keys.sort()
    for key in Keys:
        if("n" in key):
            del weightDict[key]
        elif(key in weightDict):
            revComp = gkmer_methods.getComplement(key)
            if((revComp in weightDict) and (not key == revComp)):
                weightDict[key] = weightDict[key] + weightDict[revComp]
                del weightDict[revComp]
    print("Number of unique words in dictionary after removing reverse compliments: " + str(len(weightDict)))
	return weightDict

'''
Given the dictionary where the gapped k-mer is the key, and the weight is the value,
filter out all unimportant gapped k-mers.

Inputs
weightDict:		Dictionary that contains all of the gapped k-mers as keys, and their
				weights as values.
zScore:			Threshold zScore. Anything below this threshold will not be kept.
positiveFlag:	Boolean. If TRUE, only positive zScores above the threshold are kept.
				Otherwise, all gapped k-mers with absolute value zScores above the
				threshold are kept.
				
Outputs
importantWords:	Dictionary that contains the gapped k-mers that were kept after filtering
				as keys, and an index as their value. These indices will denote the 
				gapped k-mers position in the feature vector.
'''
def filterTopGkmers(weightDict, zScore, positiveFlag):
	print("Filtering Top Gapped k-mers.")
	if positiveFlag:
		print("Flag has been set to only consider positive weights.")
		
    Values = np.asarray(weightDict.values(), dtype = 'float')
    std=np.std(Values)
    mean=np.mean(Values)
    Keys = weightDict.keys()
    Keys.sort()
    importantWords = {}
    index = 0
    for key in Keys:
        z =(weightDict.get(key)-mean)/std
        if(not positiveFlag):
            z = abs(z)
        if z > zScore:
            importantWords[key] = index
            index = index+1
    print("Number of important words: " + str(len(importantWords)))
    return importantWords
    