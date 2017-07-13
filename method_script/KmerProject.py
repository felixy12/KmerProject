import numpy as np
import math
import gkmer_methods
import write_methods
def main():
    print('Step 1: Generating gapped k-mer weights.')
    weightDict = gkmer_methods.gkmerWeight(10,4,'gm12878Data_new/gm12878_shared_ns1_svseq.fa', 'gm12878Data_new/gm12878_shared_ns1_svalpha.out', True)
    importantWords = gkmer_methods.filterTopGkmers(weightDict, 2, True)
    write_methods.writeGkmerWeights('gm12878Data_new/gm12878_top_gkmer_weights2.out', importantWords)
    print('Step 2: Creating feature matrix.')
    importantWordsIndex = gkmer_methods.convertWeightToIndex(importantWords)
    seqNameList, ftMat = gkmer_methods.getFeatures(10, 4, 'gm12878Data_new/gm12878_shared.fa', importantWordsIndex, True)
    write_methods.writeFtMat('gm12878Data_new/gm12878_shared_feature_matrix.out', ftMat, seqNameList)
    print('Step 3: Generating l-mer weights.')
    weightsList = gkmer_methods.generateLmerWeights(10, 4, importantWords, True)
    write_methods.writeLmerWeights('gm12878Data_new/gm12878_lmer_weights2.out', weightsList)
    print('Step 4: Creating position weights.')
    weightsDict = gkmer_methods.convertListToDict(weightsList)
    seqNameList, PW = gkmer_methods.positionWeights(10, weightsDict, 'gm12878Data_new/gm12878_shared.fa')
    write_methods.writePW('gm12878Data_new/gm12878_position_weights2.out', PW, seqNameList)

main()
