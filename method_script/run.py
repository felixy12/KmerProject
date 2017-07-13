def main():
    subLength = int(sys.argv[1])
    numGap = int(sys.argv[1]) - int(sys.argv[2])
    thresh = float(sys.argv[3])
    prefix = sys.argv[4]
    directoryPath = sys.argv[5]+"/"
    otherSeqFileName = sys.argv[6]
    flag = eval(sys.argv[7])

    svSeqFileName = directoryPath+prefix+"svseq.fa"
    svAlphaFileName = directoryPath+prefix+"svalpha.out"
    otherSeqFileName = directoryPath+otherSeqFileName
    if flag == 1:
        gkmerWeightsOut = directoryPath+prefix+"top_gkmer_weights"+str(thresh)+"pos.out"
    else:
        gkmerWeightsOut = directoryPath+prefix+"top_gkmer_weights"+str(thresh)+".out"
    featMatOut = directoryPath+prefix+"feature_matrix"+str(thresh)+".out"
    lmerWeightsOut = directoryPath+prefix+"lmer_weights"+str(thresh)+".out"
    posWeightsOut = directoryPath+prefix+"position_weights"+str(thresh)+".out"

    print('Step 1: Generating gapped k-mer weights.')
    weightDict = gkmer_methods.gkmerWeight(subLength, numGap, svSeqFileName, svAlphaFileName, True)
    importantWords = gkmer_methods.filterTopGkmers(weightDict, thresh, True)
    write_methods.writeGkmerWeights(gkmerWeightsOut, importantWords)
    print('Step 2: Creating feature matrix.')
    importantWordsIndex = gkmer_methods.convertWeightToIndex(importantWords)
    seqNameList, ftMat = gkmer_methods.getFeatures(subLength, numGap, otherSeqFileName, importantWordsIndex, True)
    write_methods.writeFtMat(featMatOut, ftMat, seqNameList)
    print('Step 3: Generating l-mer weights.')
    weightsList = gkmer_methods.generateLmerWeights(subLength, numGap, importantWords, True)
    write_methods.writeLmerWeights(lmerWeightsOut, weightsList)
    print('Step 4: Creating position weights.')
    weightsDict = gkmer_methods.convertListToDict(weightsList)
    seqNameList, PW = gkmer_methods.positionWeights(subLength, weightsDict, otherSeqFileName)
    write_methods.writePW(posWeightsOut, PW, seqNameList)



    

main()
