#Programmer: Chris Tralie
#Purpose: To generate series of sparse complexes for different test point cloud and
#to call Rann's pipeline to compute persistence diagrams.  Report diagrams, timings,
#and interleaving distances, as computed by wrapping around Dionysus's bottleneck
#distance binary
import subprocess
import os
import numpy as np
import time
from TestPointClouds import *
from SparseEdgeList import *

def plotDGM(dgm, color = 'b', sz = 20, label = 'dgm'):
    if dgm.size == 0:
        return
    # Create Lists
    X = list(zip(*dgm)[0]);
    Y = list(zip(*dgm)[1]);
    # set axis values
    axMin = min(min(X),min(Y));
    axMax = max(max(X),max(Y));
    axRange = axMax-axMin;
    # plot points
    plt.scatter(X, Y, sz, color,label=label)
    plt.hold(True)
    # plot line
    plt.plot([axMin-axRange/5,axMax+axRange/5], [axMin-axRange/5, axMax+axRange/5],'k');
    # adjust axis
    #plt.axis([axMin-axRange/5,axMax+axRange/5, axMin-axRange/5, axMax+axRange/5])
    # add labels
    plt.xlabel('Time of Birth')
    plt.ylabel('Time of Death')

def plot2DGMs(P1, P2, l1 = 'Diagram 1', l2 = 'Diagram 2'):
    plotDGM(P1, 'r', 10, label = l1)
    plt.hold(True)
    plt.plot(P2[:, 0], P2[:, 1], 'bx', label = l2)
    plt.legend()
    plt.xlabel("Birth Time")
    plt.ylabel("Death Time")

def savePD(filename, I):
    if os.path.exists(filename):
        os.remove(filename)
    fout = open(filename, "w")
    for i in range(I.shape[0]):
        fout.write("%g %g"%(I[i, 0], I[i, 1]))
        if i < I.shape[0]-1:
            fout.write("\n")
    fout.close()

#Wrap around Dionysus's bottleneck distance after taking the log
def getInterleavingDist(PD1, PD2):
    savePD("PD1.txt", np.log(PD1))
    savePD("PD2.txt", np.log(PD2))
    proc = subprocess.Popen(["./bottleneck", "PD1.txt", "PD2.txt"], stdout=subprocess.PIPE)
    lnd = float(proc.stdout.readline())
    return np.exp(lnd) - 1.0 #Interleaving dist is 1 + eps

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

def doRandCircleTest():
    np.random.seed(100)
    X = getRandCircle(100)
    N = X.shape[0]
    (I, J, D) = makeComplex(X, 0)
    PDsFull = getPDs(I, J, D, N, 3)
    
    (I, J, D) = makeComplex(X, 0.1)
    PDsApprox = getPDs(I, J, D, N, 3)
    
    print "Interleaving Dist Eps: ", getInterleavingDist(PDsFull[1], PDsApprox[1])
    print "Bottleneck Dist: ", getBottleneckDist(PDsFull[1], PDsApprox[1])
    [b, d, b1, d1] = [PDsFull[1][0, 0], PDsFull[1][0, 1], I[0, 0], I[0, 1]]
    print "b/b1 = %g/%g = "%(b, b1), b/b1
    print "b1/b = %g/%g = "%(b1, b), b1/b
    print "d/d1 = %g/%g = "%(d, d1), d/d1
    print "d1/d = %g/%g = "%(d1, d), d1/d
    
    plt.subplot(121)
    plt.plot(X[:, 0], X[:, 1], '.')
    plt.title('Random Point Cloud')
    plt.subplot(122)
    plotDGM(PDsFull[1], 'r', 10)
    I = PDsApprox[1]
    plt.plot(I[:, 0], I[:, 1], 'bx')
    plt.title('Persistence Diagram')
    plt.show()

def doBatchTestsShape(X, dim, TestName, doIntDist = True):
    N = X.shape[0]
    eps = np.linspace(0, 0.9, 46)
    NEps = len(eps)
    nedges = np.zeros(NEps)
    times = np.zeros(NEps)
    intdists = np.zeros(NEps)
    AllPDs = []
    for i in range(len(eps)):    
        (I, J, D) = makeComplex(X, eps[i])
        nedges[i] = len(I)
        tic = time.time()
        PDs = getPDs(I, J, D, N, dim+2)
        toc = time.time()
        times[i] = toc - tic
        AllPDs.append(PDs[dim])
        if i > 0 and doIntDist:
            intdists[i] = getInterleavingDist(AllPDs[0], AllPDs[i])
        plt.clf()
        plot2DGMs(AllPDs[0], AllPDs[i], 'Original', 'eps=%g'%eps[i])
        plt.title("%g%% Edges, dist = %g"%(float(nedges[i])/nedges[0]*100.0, intdists[i]))
        plt.savefig("%s%i.png"%(TestName, i), dpi=150, bbox_inches='tight')
    #Make video
    subprocess.call(["avconv", "-r", "2", "-i", "%s%%%s.png"%(TestName, "d"), "-r", "2", "-b", "30000k", "%s.ogg"%TestName])
    
    
    #Now plot the results
    plt.figure(figsize=(16, 16))
    plt.subplot(221)
    plt.plot(nedges, times)
    plt.xlabel("Number of Edges")
    plt.ylabel("Time (Seconds)")
    plt.title("Runtimes")
    
    plt.subplot(222)
    plt.plot(nedges, eps, 'b', label = 'eps')
    plt.hold(True)
    plt.plot(nedges, intdists, 'r', label = 'intdist')
    plt.xlabel('Number of Edges')
    plt.ylabel('Epsilon, Interleaving Dist')
    plt.legend()
    plt.title("Approximation Quality")
    
    ax = plt.subplot(223)
    rg = [X.min() - 0.1, X.max() + 0.1]
    if X.shape[1] == 2:
        plt.plot(X[:, 0], X[:, 1], 'b.')
        ax.set_xlim(rg[0], rg[1])
        ax.set_ylim(rg[0], rg[1])
    else:
        mF = plt.gca(projection='3d')
        mF.plot(X[:, 0], X[:, 1], X[:, 2], 'b.')
        mF.set_xlim(rg[0], rg[1])
        mF.set_ylim(rg[0], rg[1])
        mF.set_zlim(rg[0], rg[1])
    plt.title('Point Cloud')
    
    plt.subplot(224)
    i = 10
    plot2DGMs(AllPDs[0], AllPDs[i], 'Original', 'eps=%g'%eps[i])
    plt.title("Example, %g%% Edges, IntDist = %g"%(float(nedges[i])/nedges[0]*100.0, intdists[i]))
    plt.savefig("%s.png"%TestName, dpi=150, bbox_inches='tight')

    
if __name__ == '__main__':
    np.random.seed(100) #For repeatability

    X = getRandCircle(200)
    X = X + 0.1*np.random.randn(200, 1)
    doBatchTestsShape(X, 1, "Circle200Noisy")
    
    X = getFigure8(200)
    X = X + 0.1*np.random.randn(200, 1)
    doBatchTestsShape(X, 1, "Figure8Noisy")
    
    X = getRandSphere(200)
    doBatchTestsShape(X, 2, "Sphere200")
    
    X = getRandTorus(300, 4, 2)
    doBatchTestsShape(X, 1, "Torus300_H1")
    doBatchTestsShape(X, 2, "Torus300_H2")
        
