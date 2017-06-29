import itertools
import sys
import numpy as np
import math
import time

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
Loads in a dictionary with gapped k-mers or l-mers as keys, and weights 
(or index position) as the values.

Inputs
fileName:	Path that leads to the file with the words and weights.
useIndices:	Boolean. If TRUE, index position will be loaded in as values instead
			of the weights. Indices will be assigned based on order in which they
			appear in the file.
'''
def loadWeights(fileName, useIndices):
    weightsFile = open(fileName)
    kmerDict = {}
    index = 0
    for line in weightsFile:
        splitLine = line.strip().split()
        kmer = splitLine[0]
        weight = eval(splitLine[1])
        if useIndices:
        	kmerDict[kmer] = index
    	else:
        	kmerDict[kmer] = weight
        index += 1
    weightsFile.close()
    return kmerDict    
    
'''
Creates the reverse complement of a given DNA sequence. If a gapped sequence is given,
the gaps are ignored (they stay as gaps).

Inputs
seq:		A string containing the DNA sequence

Outputs
seq:		A string containing the reverse complement
'''
def getComplement(seq):
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
            countDict[gkmer]=countDict.get(gkmer,0.0)+1
            
    Values = np.square(countDict.values())    
    normFactor = np.sqrt(np.sum(Values))
    for key in countDict.keys():
        countDict[key] = (countDict[key]/normFactor)*alpha
    return countDict
    

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
    seqNameList, seqList = loadFasta(seqFilePath)
    alphaList = loadAlpha(alphaFilePath)
    print("Found " + str(len(seqNameList)) + " sequences to process.")
    print("Generating word weights.")
    for i in range(len(alphaList)):
        currSeq = i+1
        alpha = alphaList[i]
        seqLine = seqList[i]
        tempDict = gkmerWeightFromSeq(seqLine, subLength, numGap, alpha)
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
            revComp = getComplement(key)
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
def featureVector(seq, subLength, numGap, importantWords):
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

'''
Uses the dictionary of important words, along with a new fasta file in order to 
create a matrix where the rows denote sequences, and columns denote the feature counts.

Inputs
subLength:		Length of the gapped k-mer, l
numGap:			Number of gaps in the gapped k-mer. k = subLength - numGap
seqFilePath:	Path to the fasta file that contains the sequences we are generating
				feature vectors for.
importantWords:	Dictionary that contains gapped k-mers that exceed a certain zScore, and
				their index within the feature vector.
verbose:		Boolean flag that determins whether certain strings are printed.

Outputs
ftMat:			Feature matrix. Each row denotes a new sequence, and the column values
				represent the number of times that specific gapped k-mer appears in the
				sequence. The order for the gapped k-mers matches the index given in
				importantWords.
'''
def getFeatures(subLength, numGap, seqFilePath, importantWords, verbose):
    seqNameList, seqList = loadFasta(seqFilePath)

    start = time.time()
    print("Creating feature matrix.")
    ftMat = np.zeros((len(seqList), len(importantWords)), dtype = 'int8')
    for seqInd in range(len(seqList)):
    	curSeq = seqInd+1
        ftMat[seqInd,:] = featureVector(seqList[seqInd], subLength, numGap, importantWords)
        if(curSeq%500 == 0 and verbose):
            print("Currently processing sequence: " + str(seqInd))
    end = time.time()
    print("Done. Time Elapsed: " + str(int((end-start)/60)) + " minutes " + str(int((end-start)%60)) + " seconds") 
    return ftMat
   
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
    
    print("Calculating l-mer weights.")
    weightList = list()
    i = 1
    for lmer in lmerList:
        if(i%50000 == 0 and verbose):
            print("Currently on lmer: " + str(i))
        weight = 0
        kmerComp = getComplement(lmer)
        gkmerSet = set(gkmers(lmer, lmerLength, numGaps))
        gkmerSet.update(gkmers(lmerComp, lmerLength, numGaps))
        for gkmer in gkmerSet:
            if(gkmer in gkmerWeightsDict):
                weight += gkmerWeightsDict[gkmer]
        weightList.append((lmer,weight))
        i += 1
        
    weightList = sorted(weightList,key = lambda x:x[1],reverse=True)
    return weightList

'''
Uses existing defined l-mer weights to assign every l-mer in a fasta file a weight.

Inputs
subLength:		the length of the subsequence, or l.
lmerWeightDict:	Dictionary with the l-mer as the key, and the weight as the value.
fastaFilePath:	Path to the file that contains the sequences we want to generate the 
				weights for.
				
Output
PW:				Unlabeled Matrix that contains the sequence on the rows, and the weight at 
				each position in the column. There are n rows, where n is the number of
				sequences, and g columns, where g is the length of the longest sequence
				in the file subtracted by l.
'''
def positionWeights(subLength, lmerWeightDict, fastaFilePath):
    print('Loading in sequences.')
    seqNameList, seqList = loadFasta(fastaFilePath)
    longestSeq = 0
    for seq in seqList:
        if(longestSeq < len(seq)):
            longestSeq = len(seq)
    print('Number of sequences loaded in: ' + str(len(seqList)) + '\nLength of longest sequence: ' + str(longestSeq))
    
    print('Creating Position Weights.') 
    totalStart = time.time()
    start = time.time()
    PW = np.zeros((len(seqList), longestSeq-subLength+1))
    for i in range(len(seqList)):
        if((i+1) % 1000 == 0):
            end = time.time()
            print('Currently process sequence: ' + str(i+1)+"\t Time Elapsed: " + str((end-start)-(end-start)%0.01) + " seconds")
            start = end
        seq = seqList[i]
        for pos in range(len(seq)-subLength+1):
            subseq = seq[pos:pos+subLength]
            if subseq in lmerWeightDict:
        	    PW[i,pos] = lmerWeightDict[subseq]
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(seqList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds")
    return PW

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
    seqNameList, seqList = loadFasta(fastaFilePath)
    print('Loading in l-mer weights.')
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
   
