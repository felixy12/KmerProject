import itertools
import sys
import numpy as np
import math

#Given a string that points to the fasta file, output two lists, one with
#The sequence names, and the other with the actual sequences.
def loadFasta(fileName):
    seqFile = open(fileName, 'r')
    seqNameList = list()
    seqList = list()
    seq = seqFile.readline()
    while(">" not in seq and not seq==""):
        seq = seqFile.readline()
    seqNameList.append(seq.strip()[1:])
    seqLine = ""
    for seq in seqFile:
        if(">" not in seq and not seq==""): 
            seqLine = seqLine+seq.strip()
        elif(">" in seq):
            seqNameList.append(seq.strip()[1:])
            seqLine = seqLine.upper()
            seqList.append(seqLine)
            seqLine = ""
    seqLine = seqLine.upper()
    seqList.append(seqLine)
    seqFile.close()
    return seqNameList, seqList

def loadAlpha(fileName):
    alphaFile = open(fileName, 'r')
    alphaList = list()
    for aLine in alphaFile:
        aLine = aLine.strip().split()
        alpha = float(aLine[1])
        alphaList.append(alpha)
    alphaFile.close()
    return alphaList

#creates list of subStrings given the original string and the 
#desired length of each substring
def createSub(seq, subLength):
    seqList = list()
    for i in range(0, len(seq)-subLength+1):
        seqList.append(seq[i: i+subLength])
    return seqList

#Given a list of indices and the sequence, replace all char
#at the indices in seq with '-'
def replaceAtIndices(seq, indList):
    for i in indList:
        seq = seq[:i] + '-' + seq[i + 1:]
    return seq

#The desired function, uses itertools.combinations to get
#list of desired indices.
def gkmers(seq, subLength, numGap):
    subList = createSub(seq, subLength)
    comb = itertools.combinations(range(0,subLength), numGap)
    indList = list(comb)
    gappedList = list()
    for subS in subList:
        for iComb in indList:
            gappedList.append(replaceAtIndices(subS, iComb))
    return gappedList

#Quick method to get normalized kmerCounts, given the alpha weight of the sequence as well 
def getCounts(kmerList, alpha):
    countDict = {}
    for kmer in kmerList:
        countDict[kmer]=countDict.get(kmer, 0.0) + 1.0
    
    Values = np.square(countDict.values())
    normFactor = np.sqrt(np.sum(Values))
    for key in countDict.keys():
        countDict[key] = (countDict[key]/normFactor)*alpha
    return countDict

#Gets reverse compliment string
def getComplement(seq):
    seq = seq[::-1]
    seq = seq.replace('A', 'a')
    seq = seq.replace('T', 'A')
    seq = seq.replace('a', 'T')
    seq = seq.replace('C', 'c')
    seq = seq.replace('G', 'C')
    seq = seq.replace('c', 'G')
    return seq

#Gets kmerCounts on important kmers only
def getImportantCounts(kmerList, importantWords):
    counts = np.zeros(len(importantWords), dtype = 'int8')
    for kmer in kmerList:
        if(kmer in importantWords):
            index = importantWords.get(kmer)
            counts[index] = counts[index]+1
    return counts

#Program that can be run as a test
def main():
    seqNameList, seqList = loadFasta('data/ctcfData/ctcf_1x_svseq.fa')
    alphaList = loadAlpha('data/ctcfData/ctcf_1x_svalpha.out')
    print(len(seqNameList))
    print(len(seqList))
    print(len(alphaList))
