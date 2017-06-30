import numpy as np
import math
import gkmer_methods
import write_methods
def main():
    print('Step 1: Generating gapped k-mer weights.')
    weightDict = gkmer_methods.gkmerWeight(10,4,'../testData/test_svseq.fa', '../testData/test_svalpha.out', True)
    importantWords = gkmer_methods.filterTopGkmers(weightDict, 2, True)
    write_methods.writeGkmerWeights('../testData/test_top_gkmer_weights2.out', importantWords, weightDict)
    print('Step 2: Creating feature matrix.')
    importantWordsIndex = gkmer_methods.convertWeightToIndex(importantWords)
    seqNameList, ftMat = gkmer_methods.getFeatures(10, 4, '../testData/test_shared.fa', importantWordsIndex, True)
    write_methods.writeFtMat('../testData/test_shared_feature_matrix.out', ftMat, seqNameList)
    print('Step 3: Generating l-mer weights.')
    weightsList = gkmer_methods.generateLmerWeights(10, 4, importantWords, True)
    write_methods.writeLmerWeights('../testData/lmer_weights2.out', weightsList)
    print('Step 4: Creating position weights.')
    weightDict = gkmer_methods.convertListToDict(weightsList)
    seqNameList, PW = gkmer_methods.positionWeights(10, weightsDict, '../testData/test_shared.fa')
    write_methods.writePW('../testData/position_weights2.out', PW, seqNameList)

main()
