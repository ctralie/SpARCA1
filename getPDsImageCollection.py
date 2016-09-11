#Programmer: Chris Tralie
#Purpose: To compute persistent homology on an image collection, considering the
#Euclidean distance in pixel space
import numpy as np
import scipy.misc
import scipy.io as sio
import matplotlib.pyplot as plt
from TDA import *
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: python getPDsImageCollection.py <directory name>"
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
    
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0 #Numerical precision
    D = np.sqrt(D)
    sio.savemat("D.mat", {"D":D})
    
    
    #Step 2: Compute persistence diagrams
    print "Doing filtration..."
    PDs = doRipsFiltrationDM(D, 2)
    plotDGM(PDs[1], color = 'b')
    plt.hold(True)
    plotDGM(PDs[2], color = 'r')
#    plotDGM(PDs[3], color = 'k')
    plt.title('Persistence Diagrams')
    plt.show()
