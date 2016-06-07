import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def getRandCircle(N):
    X = np.random.randn(N, 2)
    X = X/np.sqrt(np.sum(X**2, 1))[:, None]
    return X

def getRandSphere(N):
    X = np.random.randn(N, 3)
    X = X/np.sqrt(np.sum(X**2, 1))[:, None]
    return X

def getRandTorus(N, R, r):
    theta = 2*np.pi*np.random.rand(N)
    phi = 2*np.pi*np.random.rand(N)
    X = np.zeros((N, 3))
    X[:, 0] = (R + r*np.cos(theta))*np.cos(phi)
    X[:, 1] = (R + r*np.cos(theta))*np.sin(phi)
    X[:, 2] = r*np.sin(theta)
    return X

if __name__ == '__main__':
    X = getRandTorus(200, 2, 0.5)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(X[:, 0], X[:, 1], X[:, 2], '.')
    plt.show()
