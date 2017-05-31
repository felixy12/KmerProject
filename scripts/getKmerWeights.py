#Script that generates the weight for each 10mer.
import operator
import gkmerHelpers
import sys

def generateKmer(kmerList, k):
    NTList = ['T','C','G','A']
    if(k==0):
        return kmerList
    kmerListNew = list()
    for kmer in kmerList:
        for NT in NTList:
            kmerListNew.append(kmer+NT)
    return generateKmer(kmerListNew, k-1)



def main():
    gkmerLength = eval(sys.argv[1])
    numGaps = eval(sys.argv[1])-eval(sys.argv[2])
    directory = sys.argv[3]
    weightsFileName = sys.argv[4]
    outputFileName = sys.argv[5]
    kmerList = [""]

    print("Loading in gkmer weights")
    weightsFile = open(directory + '/' + weightsFileName)
    gkmerWeightsDict = {}
    for line in weightsFile:
        splitLine = line.strip().split()
        gkmer = splitLine[0]
        weight = eval(splitLine[1])
        gkmerWeightsDict[gkmer] = weight
    weightsFile.close()
    print('Number of gkmers loaded in: ' + str(len(gkmerWeightsDict)))
    
    kmerList = [""]
    kmerList = generateKmer(kmerList,10)
    
    print("Generating kmer weights.")
    weightList = list()
    i = 1
    for kmer in kmerList:
        if(i%10000 == 0):
            print("Currently on kmer: " + str(i))
        weight = 0
        kmerComp = gkmerHelpers.getComplement(kmer)
        gkmerSet = set(gkmerHelpers.gkmers(kmer, gkmerLength, numGaps))
        gkmerSet.update(gkmerHelpers.gkmers(kmerComp, gkmerLength, numGaps))
        for gkmer in gkmerSet:
            if(gkmer in gkmerWeightsDict):
                weight += gkmerWeightsDict[gkmer]
        weightList.append((kmer,weight))
        i += 1
    print(len(weightList))  
    weightList = sorted(weightList,key = lambda x:x[1],reverse=True)
    outFile = open(directory+'/'+outputFileName,'w')
    print("Writing results to file.")
    for tup in weightList:
        outFile.write(tup[0]+'\t'+str(tup[1])+'\n')
    outFile.close()
main()
