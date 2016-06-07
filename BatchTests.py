import subprocess
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

#Wrap around Dionysus's bottleneck distance after taking the log
def getInterleavingDist(PD1, PD2):    
    np.savetxt("PD1.txt", np.log(PD1))
    np.savetxt("PD2.txt", np.log(PD2))
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    lnd = float(x.stdout.readline())
    return np.exp(lnd)

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
    writeResults(I, J, D, N, "temp.dimacs")
    subprocess.call(["./phatclique", "-i", "temp.dimacs", "-m", "%i"%m, "-o" "temp.results"])
    return parsePDs("temp.results")

if __name__ == '__main__':
    np.random.seed(100)
    X = getRandCircle(100)
    N = X.shape[0]
    (I, J, D) = makeComplex(X, 0)
    PDs = getPDs(I, J, D, N, 4)
    
    plt.subplot(121)
    plt.plot(X[:, 0], X[:, 1], '.')
    plt.title('Random Point Cloud')
    plt.subplot(122)
    plotDGM(PDs[1])
    plt.title('Persistence Diagram')
    plt.show()
