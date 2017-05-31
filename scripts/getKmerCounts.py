import numpy
import sys
import gkmerHelpers
def main():
    gkmerLength = eval(sys.argv[1])
    numGaps = eval(sys.argv[1])-eval(sys.argv[2])
    directory = sys.argv[3]
    fastaFileName = sys.argv[4]
    kmerFileName = sys.argv[5]
    outputFileName = sys.argv[6]

    print('Loading in sequences.')
    seqNameList, seqList = gkmerHelpers.loadFasta(directory+'/'+fastaFileName)
    print('Loading in k-mer weights.')
    kmerFile = open(directory+'/'+kmerFileName)
    
    kmerCountDict = {}
    kmerOrder = list()
    for line in kmerFile:
        splitLine = line.strip().split()
        kmerOrder.append(splitLine[0])
        kmerCountDict[splitLine[0]] = 0
    for i in range(len(seqList)):
        if((i+1) % 1000 == 0): 
            print('Currently processing sequence: ' + str(i+1))
        seq = seqList[i]
        kmerList = gkmerHelpers.createSub(seq, gkmerLength)
        for k in kmerList:
            kmerCountDict[k]+=1
        i += 1
    print('Saving kmer counts to file: ' + directory+'/'+outputFileName)
    outputFile = open(directory+'/'+outputFileName,'w')
    for kmer in kmerOrder:
        outputFile.write(kmer+'\t'+str(kmerCountDict[kmer])+'\n')
    outputFile.close()
main()
