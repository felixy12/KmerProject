import numpy as np
#import matplotlib.pyplot as plt
#import pylab

def calculateThresh(posArray, negArray):
    posArray = np.sort(posArray)
    negArray = np.sort(negArray)
    index = len(negArray)-1
    stop = False
    while not stop:
        thresh = negArray[index]
        numPos = np.sum(posArray>=thresh)
        numNeg = np.sum(negArray>=thresh)
        print("NumPos: "+str(numPos)+"\tNumNeg: "+str(numNeg))
        index=index-1
        #print("index: "+str(index)+"\tThreshold: "+str(thresh))
        #if(numPos/numNeg > 100):
        #    return 

def loadData(fileName):
    data = np.loadtxt(fileName, dtype='S')
    data = np.delete(data,[0,291,292,293,294,295,296,297,298,299],axis=1)
    data = data.astype('f')
    return data

def writeVector(weightVector, fileName):
    out = open(fileName, 'w')
    toWrite = ""
    for i in range(len(weightVector)):
        toWrite += str(weightVector[i]) + "\t"
    toWrite = toWrite.strip()+'\n'
    out.write(toWrite)
    out.close()            

def main():
    print('Reading in positive matrix')
    pos_top_file = 'pos_top_sum_positive.out'
    PTM = loadData(pos_top_file)
    PTV = PTM.flatten()
    print('Reading in negative matrix')
    neg_top_file = 'pos_top_sum_negative.out'
    NTM = loadData(neg_top_file)
    NTV = NTM.flatten()
    print(len(PTV))
    print(len(NTV))
    #writeVector(PTV, 'positive_top_vector.out')
    #writeVector(NTV, 'negative_top_vector.out')
    plt.figure(1)
    plt.subplot(211)
    plt.hist(PTV,range=(0,100))
    plt.subplot(212)
    plt.hist(NTV,range=(0,100))
    pylab.show()
    #calculateThresh(PTV,NTV)
    #writeVector(PTV, 'top_positive_vector.out')

main()
