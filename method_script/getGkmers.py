import gkmerHelpers
import sys 
import time
import numpy as np
#To run, use >python getGkmer.py <length of subsequence> <number of non-gaps> <z-score threshold> <filename prefix (in this case, ctcf_1x_)> <path to directory with data> <flag>
def gkmerWeight(subLength, numGap, zScore, seqFilePath, alphaFilePath, positiveFlag, verbose):
        
    #variables used to analyze script runtime and size constraints.
    strLength = 0;
    dictSize = 0;
    totalStart = time.time()
    start = totalStart
    end = start
    seqDict = {}
    numPerAnalysis = 500 
    numSeqs = 0    

    #Meat of the script
    #This next line assumes the first line is a sequence position, and no header lines exist. If there are header lines, please remove them.
    numSeqs = numSeqs+1
    print("Calculating gapped k-mer weights")
    if(positiveFlag == 1):
        print("Flag has been set to only consider positive top weights.")
    print("Loading files.")
    seqNameList, seqList = gkmerHelpers.loadFasta(seqFilePath)
    alphaList = gkmerHelpers.loadAlpha(alphaFilePath)
    print("Found " + str(len(seqNameList)) + " sequences to process.")
    print("Generating word weights.")
    for i in range(len(alphaList)):
        currSeq = i+1
        alpha = alphaList[i]
        seqLine = seqList[i]
        #Generating results and analyzing runtime.
        words = gkmerHelpers.gkmers(seqLine, subLength, numGap)
        tempDict = gkmerHelpers.getCounts(words, alpha)
        dictSize = dictSize + len(tempDict)
        for key in tempDict.keys():
            seqDict[key] = seqDict.get(key, 0) + tempDict[key]
        strLength = strLength+len(seqLine)
        if(currSeq%numPerAnalysis == 0 and currSeq != 0): 
            end = time.time()
            print("Summary of sequences " + str(currSeq-numPerAnalysis+1) + "-"+str(currSeq))    
            print("Time Elapsed: " + str(end-start) + " seconds")
            print("Average length of sequence: " + str(strLength/numPerAnalysis)+"bp")
            print("Average size of dictionary for each sequence: " + str(dictSize/numPerAnalysis) + " entries")
            print("Current number of unique words: " + str(len(seqDict.keys())))
            print("")
            strLength = 0
            dictSize = 0;
            start = end
    totalEnd = time.time()
    print("Done generating. Total time elapsed to process " + str(len(alphaList)) + " sequences: " + str(int((totalEnd-totalStart)/60)) + " minutes " + str(int((totalEnd-totalStart)%60)) + " seconds")

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
        z =(seqDict.get(key)-mean)/std
        if(not flag == 1):
            z = abs(z)
        if z > zScore:
            importantWords[key] = index
            index = index+1
    print("Number of important words: " + str(len(importantWords)))
    outFileName = directoryPath+prefix+"topweights"+str(zScore)
    if(flag == 1):
        outFileName = outFileName+"_pos"
    wordsFile = open(outFileName+".out", 'w')

    Keys = importantWords.keys()
    Keys.sort()
    for word in Keys:
        wordsFile.write(word+"\t"+str(seqDict.get(word))+"\n")
    wordsFile.close()

main()
