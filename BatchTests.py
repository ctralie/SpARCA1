import subprocess
import os
import numpy as np
from TestPointClouds import *
from SparseEdgeList import *

def plotDGM(dgm, color = 'b', sz = 20):
    # Create Lists
    X = list(zip(*dgm)[0]);
    Y = list(zip(*dgm)[1]);
    # set axis values
    axMin = min(min(X),min(Y));
    axMax = max(max(X),max(Y));
    axRange = axMax-axMin;
    # plot points
    plt.scatter(X, Y, sz, color)
    plt.hold(True)
    # plot line
    plt.plot([axMin-axRange/5,axMax+axRange/5], [axMin-axRange/5, axMax+axRange/5],'k');
    # adjust axis
    plt.axis([axMin-axRange/5,axMax+axRange/5, axMin-axRange/5, axMax+axRange/5])
    # add labels
    plt.xlabel('Time of Birth')
    plt.ylabel('Time of Death')

def savePD(filename, I):
    fout = open(filename, "w")
    for i in range(I.shape[0]):
        fout.write("%g %g"%(X[i, 0], X[i, 1]))
        if i < I.shape[0]-1:
            fout.write("\n")
    fout.close()

#Wrap around Dionysus's bottleneck distance after taking the log
def getInterleavingDist(PD1, PD2):
    savePD("PD1.txt", np.log(PD1))
    savePD("PD2.txt", np.log(PD2))
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    lnd = float(proc.stdout.readline())
    return np.exp(lnd)

def getBottleneckDist(PD1, PD2):
    savePD("PD1.txt", PD1)
    savePD("PD2.txt", PD2)
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    return float(proc.stdout.readline())

def parsePDs(filename):
    PDs = {}
    fin = open(filename)
    for l in fin.readlines():
        fs = [float(s.rstrip()) for s in l.split()]
        dim = int(fs[0])
        if not dim in PDs:
            PDs[dim] = []
        if fs[-2] == fs[-1]:
            continue #Don't display classes which die instantly
        PDs[dim].append(fs[-2:])
    fin.close()
    ret = []
    count = 0
    for i in range(len(PDs)):
        ret.append(np.array(PDs[i]))
    return ret

def getPDs(I, J, D, N, m):
    if os.path.exists("temp.dimacs"):
        os.remove("temp.dimacs")
    writeResults(I, J, D, N, "temp.dimacs")
    if os.path.exists("temp.results"):
        os.remove("temp.results")
    subprocess.call(["./phatclique", "-i", "temp.dimacs", "-m", "%i"%m, "-o" "temp.results"])
    return parsePDs("temp.results")

if __name__ == '__main__':
    np.random.seed(100)
    X = getRandCircle(100)
    N = X.shape[0]
    (I, J, D) = makeComplex(X, 0)
    PDsFull = getPDs(I, J, D, N, 3)
    
    (I, J, D) = makeComplex(X, 0.2)
    PDsApprox = getPDs(I, J, D, N, 3)
    
    plt.subplot(121)
    plt.plot(X[:, 0], X[:, 1], '.')
    plt.title('Random Point Cloud')
    plt.subplot(122)
    plotDGM(PDsFull[1], 'r', 10)
    I = PDsApprox[1]
    plt.plot(I[:, 0], I[:, 1], 'bx')
    plt.title('Persistence Diagram')
    plt.show()
