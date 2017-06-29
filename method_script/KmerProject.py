import numpy as np
import math
import gkmer_methods

def main():
    weightDict = gkmer_methods.gkmerWeight(10,4,'../testData/test_svseq.fa', '../testData/test_svalpha.out', True)
    importantWords = filterTopGkmers(weightDict, 2, True)
    writeGkmerWeights('../testData/test_top_gkmer_weights2.out', importantWords, weightDict)
main()
