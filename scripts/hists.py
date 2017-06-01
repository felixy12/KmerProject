import numpy as np
import matplotlib.pyplot as plt
import pylab
import prettyplotlib as ppl
from prettyplotlib import brewer2mpl

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
    pos_top_file = '../gm12878Data/position_weights_pos.out'
    PWPM = loadData(pos_top_file)
    print('Reading in negative matrix')
    neg_top_file = '../gm12878Data/position_weights_neg.out'
    PWNM = loadData(neg_top_file)
    
    #Process that generates histograms of kmer weights
    '''
    PWPV = PWPM.flatten()
    PWNV = PWNM.flatten()
    plt.figure(1)
    plt.subplot(211)
    plt.hist(PWPV,range=(30,100))
    plt.subplot(212)
    plt.hist(PWNV,range=(30,100))
    pylab.show()
    '''
    #Process that generates histograms of sequence weights
    pos_seq_weights = np.sum(PWPM,1)
    neg_seq_weights = np.sum(PWNM,1)
    pos_hist = np.histogram(pos_seq_weights, np.linspace(-5000, 10000, 1000))[0]
    pos_hist = pos_hist.astype('f')
    print(pos_hist)
    neg_hist = np.histogram(neg_seq_weights, np.linspace(-5000, 10000, 1000))[0]
    neg_hist = neg_hist.astype('f')
    print(neg_hist)
    print(pos_hist.shape)
    overlap = np.minimum(pos_hist, neg_hist)
    print('overlap done')
    overlap = np.sum(overlap)
    fraction_overlap = overlap/len(pos_seq_weights)
    fraction_overlap = fraction_overlap - fraction_overlap%0.0001
    print(fraction_overlap)
    plt.figure(1)
    plt.hist(pos_seq_weights, bins=100, alpha = 0.5, label = 'positive sequence weights')
    plt.hist(neg_seq_weights, bins=100, alpha = 0.5, label = 'negative sequence weights')
    plt.text(3000, 500, "Fraction Overlap: "+str(fraction_overlap))
    plt.legend(loc='upper right')
    pylab.show()
    
main()
