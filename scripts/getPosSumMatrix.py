# Run with the following Command: python getPosMatrix.py <gkmerLength> <num non gapped> <fasta file name> <top weights file name> <directory> <output file name>
import numpy as np
import gkmerHelpers
import sys
import time
import random
def getPosSumVector(seq, gkmerWeightsDict, gkmerLength, numGaps, numPos):
    numGkmers = len(gkmerWeightsDict)
    posSumVector = np.zeros(numPos, dtype='float')
    subSeqsList = gkmerHelpers.createSub(seq, gkmerLength)
    for i in range(len(subSeqsList)):
        seq = subSeqsList[i]
        seqComp = gkmerHelpers.getComplement(seq)
        gkmerSet = set(gkmerHelpers.gkmers(seq, gkmerLength, numGaps))
        gkmerSet.update(gkmerHelpers.gkmers(seqComp, gkmerLength, numGaps))
        for gkmer in gkmerSet:
            if(gkmer in gkmerWeightsDict):
                posSumVector[i] += gkmerWeightsDict[gkmer]
    return posSumVector

def writeMatrix(directory, posSumMatrix, seqNameList, outFileName):
    print('Writing matrix to output file.')
    name = directory + '/' + outFileName
    outFile = open(name, 'w')
    for i in range(posSumMatrix.shape[0]):
        row = seqNameList[i]+'\t'
        for j in range(posSumMatrix.shape[1]):
            row = row+str(posSumMatrix[i,j])+'\t'
        outFile.write(row.strip()+'\n')
    print('Done writing. File can be found in ' + name)
    outFile.close()

def main():
    gkmerLength = eval(sys.argv[1])
    numGaps = eval(sys.argv[1])-eval(sys.argv[2])
    faFileName = sys.argv[3]
    weightsFileName = sys.argv[4]
    directory = sys.argv[5]
    outFileName = sys.argv[6]
    print('Loading in sequences.')
    seqNameList, seqList = gkmerHelpers.loadFasta(directory+'/'+faFileName)
    longestSeq = 0
    for seq in seqList:
        if(longestSeq < len(seq)):
            longestSeq = len(seq)
    print('Number of sequences loaded in: ' + str(len(seqList)) + '\nLength of longest sequence: ' + str(longestSeq))
    
    print('Loading in gkmer weights.')
    weightsFile = open(directory+'/'+weightsFileName)
    gkmerWeightsDict = {}
    for line in weightsFile:
        splitLine = line.strip().split()
        gkmer = splitLine[0]
        weight = eval(splitLine[1])
        gkmerWeightsDict[gkmer] = weight
    weightsFile.close()
    print('Number of gkmers loaded in: ' + str(len(gkmerWeightsDict)))
    
    print('Creating Position Sum Matrix.') 
    totalStart = time.time()
    start = time.time()
    PSM = np.zeros((len(seqList), longestSeq))
    for i in range(len(seqList)):
        if((i+1) % 1000 == 0):
            end = time.time()
            print('Currently process sequence: ' + str(i+1)+"\t Time Elapsed: " + str((end-start)-(end-start)%0.01) + " seconds")
            start = end
        seq = seqList[i]
        PSM[i,:] = getPosSumVector(seq, gkmerWeightsDict, gkmerLength, numGaps, longestSeq)
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(seqList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds")
    writeMatrix(directory, PSM, seqNameList, outFileName)
main()

        
