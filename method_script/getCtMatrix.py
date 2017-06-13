import numpy as np
import sys
import gkmerHelpers
import time

#to run, use >python <length of subsequence> <number of non-gaps> <sequence file name> <top words file name> <directory path> <output file name>
def main():
    subLength = int(sys.argv[1])
    numGap = int(sys.argv[1])-int(sys.argv[2])
    seqFileName = sys.argv[3]
    weightsFileName = sys.argv[4]
    directory = sys.argv[5]+'/'
    outFileName = sys.argv[6]
    weightsFile = open(directory+weightsFileName, 'r')

    seqList = list()
    seqNameList, seqList = gkmerHelpers.loadFasta(directory+seqFileName)
    importantWords = {}

    i = 0
    for word in weightsFile:
        wordWeight = word.strip().split()
        importantWords[wordWeight[0]] = i
        i = i+1
    weightsFile.close()

    start = time.time()
    #Creates the count matrix
    print("Creating counts matrix.")
    ctMat = np.zeros((len(seqList), len(importantWords)), dtype = 'int8')
    for seqInd in range(len(seqList)):
        words = gkmerHelpers.gkmers(seqList[seqInd], subLength, numGap)
        ctMat[seqInd,:] = gkmerHelpers.getImportantCounts(words, importantWords)
        if(seqInd%500 == 0 and seqInd != 0):
            print("Currently processing sequence: " + str(seqInd))
    end = time.time()
    print("Done. Time Elapsed: " + str(int((end-start)/60)) + " minutes " + str(int((end-start)%60)) + " seconds") 
    print("Writing counts matrix to file " + directory+outFileName)
    #Writes the count matrix to a file
    ctFile = open(directory+outFileName, 'w')
    for seqInd in range(len(seqNameList)):
        ctFile.write(seqNameList[seqInd])
        for i in range(len(importantWords)):
            ctFile.write("\t"+str(ctMat[seqInd,i]))
        ctFile.write("\n")
    ctFile.close()
    
main()
