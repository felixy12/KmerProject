# Run with the following Command: python getPosMatrix.py <gkmerLength> <num non gapped> <num sequences> <fasta file name> <top weights file name> <directory>
import numpy as np
import gkmerHelpers
import sys
import time
import random
def getPosMatrix(seq, topWeightsList, gkmerLength, numGaps):
    numPos = len(seq)-gkmerLength+1
    numGkmers = len(topWeightsList)
    posMatrix = np.zeros((numPos, numGkmers), dtype='uint8')
    indexMapping = {}

    subSeqsList = gkmerHelpers.createSub(seq, gkmerLength)
    for i in range(len(subSeqsList)):
        seq = subSeqsList[i]
        seqComp = gkmerHelpers.getComplement(seq)
        gkmerSet = set(gkmerHelpers.gkmers(seq, gkmerLength, numGaps))
        gkmerSet.update(gkmerHelpers.gkmers(seqComp, gkmerLength, numGaps))
        for j in range(len(topWeightsList)):
            if(topWeightsList[j] in gkmerSet):
                posMatrix[i, j] = 1
    return posMatrix

def writeMatrix(directory, seqName, posMatrix, topWeightsList):
    outFile = open(directory + '/' + seqName+"_posmatrix.out", 'w')
    header = ""
    for gkmer in topWeightsList:
        header = header+gkmer+'\t'
    outFile.write(header.strip()+'\n')
    for i in range(posMatrix.shape[0]):
        row = ""
        for j in range(posMatrix.shape[1]):
            row = row+str(posMatrix[i,j])+'\t'
        outFile.write(row.strip()+'\n')
    outFile.close()

def main():
    gkmerLength = eval(sys.argv[1])
    numGaps = eval(sys.argv[1])-eval(sys.argv[2])
    numSeq = eval(sys.argv[3])
    faFileName = sys.argv[4]
    weightsFileName = sys.argv[5]
    directory = sys.argv[6]
    seqNameList, seqList = gkmerHelpers.loadFasta(directory+'/'+faFileName)
    weightsFile = open(directory+'/'+weightsFileName)
    topWeightsList = list()
    indList = list()
    numSeq = 20
    for line in weightsFile:
        gkmer = line.strip().split()[0]
        topWeightsList.append(gkmer)
    weightsFile.close()

    for i in range(numSeq):
        randInd = random.randint(0, len(seqList)-1)
        if(randInd not in indList):
            indList.append(randInd)
            seq = seqList[randInd]
            seqName = seqNameList[randInd]
            print('generating position matrix for sequence: ' + seqName)
            posMatrix = getPosMatrix(seq, topWeightsList, gkmerLength, numGaps)
            print('writing result')
            writeMatrix(directory, seqName, posMatrix, topWeightsList)
main()

        
