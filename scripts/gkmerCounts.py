#This is the old script that ran everything at once. This script is now divided into two, getGkmers.py, which outputs important gkmers into topweights.out, and
#getCounts.py, which outputs positive sequence counts int ctcounts.out

import gkmerHelpers
import sys
import time
import numpy as np
#To run, use >python gkmerWeights.py <length of subsequence> <number of gaps (k)> <z-score threshold> <filename prefix (in this case, ctcf_1x_)> <path to directory with data>
def main():
    #setting up all necessary files and function parameters read in from command line.
    subLength = int(sys.argv[1])
    numGap = int(sys.argv[1]) - int(sys.argv[2])
    zScore = float(sys.argv[3])
    prefix = sys.argv[4]
    directoryPath = sys.argv[5]+"/"
    seqFile = open(directoryPath+prefix+"svseq.fa", 'r')
    alphaFile = open(directoryPath+prefix+"svalpha.out", 'r')
    
    #variables used to analyze script runtime and size constraints.
    strLength = 0;
    dictSize = 0;
    totalStart = time.time()
    start = totalStart
    end = start
    seqDict = {}
    numPerAnalysis = 500
    
    #Lists used to generate the ctmatrix. seqNameList holds the positions of each sequence, with seqList holds the actual sequence. These pairs are matched by identical indices.
    seqNameList = list()
    seqList = list()

    print("Generating gapped word weights.")
    #Meat of the script
    #This next line assumes the first line is a sequence position, and no header lines exist. If there are header lines, please remove them.
    seq = seqFile.readline()
    while(">" not in seq and not seq==""):
        seq = seqFile.readline()
    seqNameList.append(seq)
    for aLine in alphaFile:
        #Grabbing the two corresponding alpha and seq (along with the position), and appending values onto their respective lists.
        aLine=aLine.strip().split()
        alpha = float(aLine[1])
        seq=seqFile.readline()
        seqLine = ""
        while(">" not in seq and not seq==""): 
            seqLine = seqLine+seq.strip()
            seq=seqFile.readline()
            if(">" in seq):
                seqNameList.append(seq.strip())
        seqLine = seqLine.upper()
        seqList.append(seqLine)
        #Generating results and analyzing runtime.
        words = gkmerHelpers.gkmers(seqLine, subLength, numGap)
        tempDict = gkmerHelpers.getCounts(words, alpha)
        dictSize = dictSize + len(tempDict)
        for key in tempDict.keys():
            seqDict[key] = seqDict.get(key, 0) + tempDict[key]
        strLength = strLength+len(seqLine)
        if(len(seqNameList)%numPerAnalysis == 0):
            end = time.time()
            print("Summary of sequences " + str(len(seqNameList)-numPerAnalysis+1) + "-"+str(len(seqNameList)))    
            print("Time Elapsed: " + str(end-start) + " seconds")
            print("Average length of sequence: " + str(strLength/numPerAnalysis)+"bp")
            print("Average size of dictionary for each sequence: " + str(dictSize/numPerAnalysis) + " entries")
            print("Current number of unique words: " + str(len(seqDict.keys())))
            print("")
            strLength = 0
            dictSize = 0;
            start = end
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(seqNameList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds") 
    
    #Does post processing on the matrix, including combining reverse compliments
    Keys = seqDict.keys()
    if("C--CT--TGG" in seqDict and "CCA--AG--G"):
        seqDict["C--CT--TGG"] = seqDict["C--CT--TGG"] + seqDict["CCA--AG--G"]
        del seqDict["CCA--AG--G"]
    Keys.sort()
    for key in Keys:
        if("n" in key):
            del seqDict[n]
        elif(key in seqDict):
            revComp = gkmerHelpers.getComplement(key)
            if((revComp in seqDict) and (not key == revComp)):
                seqDict[key] = seqDict[key] + seqDict[revComp]
                del seqDict[revComp]
    print("Number of unique words in dictionary after removing reverse compliments: " + str(len(seqDict))) 
    
    #Generates list of words that surpass a certain z-score.
    Values = np.asarray(seqDict.values(), dtype = 'float')
    std=np.std(Values)
    mean=np.mean(Values)
    Keys = seqDict.keys()
    Keys.sort()
    importantWords = {}
    index = 0
    for key in Keys:
        z =abs((seqDict.get(key)-mean)/std)
        if z > zScore:
            importantWords[key] = index
            index = index+1
    print("Number of important words: " + str(len(importantWords)))

    #Creates the count matrix
    print("Creating counts matrix.")
    ctMat = np.zeros((len(seqList), len(importantWords)))
    for seqInd in range(len(seqList)):
        words = gkmerHelpers.gkmers(seqList[seqInd], subLength, numGap)
        ctMat[seqInd,:] = gkmerHelpers.getImportantCounts(words, importantWords)
    end = time.time()
    print("Done. Time Elapsed: " + str(int((end-start)/60)) + " minutes " + str(int((end-start)%60)) + " seconds") 
    print("Writing counts matrix to file " + directoryPath+prefix+"ctmatrix.out")
    #Writes the count matrix to a file
    ctFile = open(directoryPath+prefix+"ctmatrix.out", 'w')
    for seqInd in range(len(seqNameList)):
        ctFile.write(seqNameList[seqInd])
        for i in range(len(importantWords)):
            ctFile.write("\t"+str(ctMat[seqInd,i]))
        ctFile.write("\n")

    print("Writing important words to file " + directoryPath+prefix+"topweights.out")
    #Writes important words to a file, along with it's respective weight
    wordsFile = open(directoryPath+prefix+"topweights.out", 'w')
    Keys = importantWords.keys()
    Keys.sort()
    for word in Keys:
        wordsFile.write(word+"\t"+str(seqDict.get(word))+"\n")

    wordsFile.close()
    ctFile.close()
    seqFile.close()
    alphaFile.close()

main()

