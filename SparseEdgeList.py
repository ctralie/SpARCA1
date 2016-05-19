#Programmer: Chris Tralie
#Purpose: Return a warped edge list based on Cavanna/Jahanseir/Sheehy's approach
#to approximate rips filtrations over some metric space, using a greedy permutation.
#The goal is to return many fewer edges than all O(N^2)

from sys import argv, exit
import numpy as np


#Purpose: Naive O(N^2) algorithm to do the greedy permutation
#Inputs: D (NxN distance matrix for points)
#Returns: (permutation (N-length array of indices), 
#lambdas (N-length array of insertion radii))
def getGreedyPerm(D):
    N = D.shape[0]
    #By default, takes the first point in the list to be the
    #first point in the permutation, but could be random
    perm = np.zeros(N, dtype=np.int64)
    lambdas = np.zeros(N)
    ds = D[0, :]
    for i in range(1, N):
        idx = np.argmax(ds)
        perm[i] = idx
        lambdas[i] = ds[idx]
        ds = np.minimum(ds, D[idx, :])
    return (perm, lambdas)

#Purpose: To return the sparse edge list with the warped distances, sorted
#by weight
#Inputs: lambdaso (insertion radii for points), eps (epsilon approximation constant),
#D (NxN distance matrix, okay to modify because last time it's used)
def getSparseEdgeList(lambdaso, eps, D):
    N = D.shape[0]
    
    #Create initial sparse list candidates (Lemma 6)
    nBounds = ((eps**2+3*eps+2)/eps)*lambdaso #Search neighborhoods    
    D[D > nBounds[:, None]] = np.inf #Set all distances outside of search neighborhood to infinity
    [I, J] = np.meshgrid(np.arange(N), np.arange(N))
    I = I[D < np.inf]
    J = J[D < np.inf]
    D = D[D < np.inf]
    
    #Prune sparse list and update warped edge lengths (Algorithm 3 pg. 14)
    lambdas = np.zeros((len(I), 2))
    lambdas[:, 0] = np.minimum(lambdaso[I], lambdaso[J])
    lambdas[:, 1] = np.maximum(lambdaso[I], lambdaso[J])
    t = np.arange(len(I))
    t = t[D <= np.sum(lambdas, 1)*(1+eps)/eps]
    I = I[t]
    J = J[t]
    D = D[t]
    lambdas = lambdas[t, :]
    t = np.ones(len(I))
    t[D <= 2*lambdas[:, 0]*(1+eps)/eps] = 0
    D[t == 1] = 2.0*(D[t == 1] - lambdas[t == 1, 0]*(1+eps)/eps) #Multiply by 2 since this is Rips not Cech
    return (I, J, D)

if __name__ == '__main__':
    #Step 1: Load in point cloud
    if len(argv) < 4:
        print "Usage: WarpedRipsMetric <point cloud filename> <epsilon> <outname> <DEBUG (optional)>"
        exit(0)
    #Assumes a point clouds is in a text file with each point on its own line
    #and dimensions separated by spaces    
    #X = sio.loadmat(argv[1])['X']
    #N = X.shape[0]
    fin = open(argv[1])
    X = fin.readlines()
    for i in range(len(X)):
        X[i] = [float(x) for x in X[i].split()]
    N = len(X)
    X = np.array(X)
    eps = float(argv[2])
    
    #Step 2: Compute all pairwise distances
    XSqr = np.sum(X**2, 1)
    D = XSqr[:, None] + XSqr[None, :] - 2*X.dot(X.T)
    D[D < 0] = 0 #Numerical precision
    D = np.sqrt(D)
    
    I = []
    J = []
    #Step 2.5: If epsilon is zero, we're done (use all edges)
    if eps == 0:
        [I, J] = np.meshgrid(np.arange(N), np.arange(N))
        idx = (I > J)
        I = I[idx == 1]
        J = J[idx == 1]
        D = D[idx == 1]
    else:
        #Step 3: Compute greedy permutation
        (perm, lambdas) = getGreedyPerm(D)
        
        #Step 4: Compute the warped edge list
        #Put the insertion radii back in the order of the array
        lambdaso = np.zeros(len(lambdas))
        lambdaso[perm] = lambdas
        (I, J, D) = getSparseEdgeList(lambdaso, eps, D)
    
    #Step 5: Write to file
    #First sort edges by weight
    idx = np.argsort(D)
    [I, J, D] = [I[idx], J[idx], D[idx]]
    fout = open(argv[3], 'w')
    fout.write("p edge %i %i\n"%(N, len(I)))
    for i in range(len(I)):
        fout.write("e %i %i %g\n"%(I[i], J[i], D[i]))
    fout.close()
