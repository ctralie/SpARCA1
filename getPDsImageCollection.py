#Programmer: Chris Tralie
#Purpose: To compute persistent homology on an image collection, considering the
#Euclidean distance in pixel space
import numpy as np
import scipy.misc
import matplotlib.pyplot as plt
from SparseEdgeList import *
from BatchTests import *
import sys

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: python getPDsImageCollection.py <directory name> <num frames>"
        sys.exit(0)
    dirName = sys.argv[1]
    N = int(sys.argv[2])
    
    #Step 1: Load image collection
    I0 = scipy.misc.imread("%s/0.png"%dirName)
    I0 = I0.flatten()
    X = np.zeros((N, len(I0)))
    X[0, :] = I0
    for i in range(1, N):
        print "Reading %i"%i
        I = scipy.misc.imread("%s/%i.png"%(dirName, i))
        X[i, :] = I.flatten()
    
    #Step 2: Compute persistence diagrams
    eps = 0
    print "Doing filtration..."
    (I, J, D) = makeComplex(X, eps)
    PDs = getPDs(I, J, D, N, 4)
    plotDGM(PDs[1])
    plt.title('1D Persistence Diagram')
    plt.show()
    plotDGM(PDs[2])
    plt.title('2D Persistence Diagram')
    plt.show()
