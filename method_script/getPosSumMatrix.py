import numpy as np
import gkmer_methods
import time

'''
Uses existing defined l-mer weights to assign every l-mer in the fasta file a weight.

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
    seqNameList, seqList = gkmerHelpers.loadFasta(fastaFilePath)
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
        	subseq = seq(pos:pos+subLength)
        	if subseq in lmerWeightDict:
        		PW[i,pos] = lmerWeightDict[subseq]
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(seqList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds")
	return PW
main()

        
