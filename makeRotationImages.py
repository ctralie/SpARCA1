#Programmer: Chris Tralie
#Purpose: To view a mesh from a bunch of random directions on the
#sphere
import sys
sys.path.append("../")
sys.path.append("../S3DGLPy")
import numpy as np
import scipy
from Primitives3D import *
from PolyMesh import *
from MeshCanvas import *
import matplotlib.pyplot as plt

#########################################################
##                UTILITY FUNCTIONS                    ##
#########################################################

class ShapeGLCanvas(BasicMeshCanvas):
    def __init__(self, parent, angles, dirName, mesh):
        BasicMeshCanvas.__init__(self, parent)
        self.angles = angles #Angles to rotate the camera through
        self.dirName = dirName
        self.frameNum = 0
        self.mesh = mesh
        self.initMeshBBox()
        self.displayMeshVertices = False
        
        self.Refresh()
    
    def setupPerspectiveMatrix(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(180.0*self.camera.yfov/M_PI, float(self.size.x)/self.size.y, 0.001, 100)
    
    def repaint(self):
        self.setupPerspectiveMatrix()
        glClearColor(0.0, 0.0, 0.0, 0.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        glEnable(GL_LIGHTING)
        
        #Rotate camera
        self.camera.theta = self.angles[self.frameNum, 0]
        self.camera.phi = self.angles[self.frameNum, 1]
        self.camera.updateVecsFromPolar()

        #Set up modelview matrix
        self.drawMeshStandard()       

        saveImageGL(self, "%s/%i.png"%(self.dirName, self.frameNum))
        if self.frameNum < self.angles.shape[0] - 1:
            self.frameNum += 1
            self.Refresh()
        else:
            self.parent.Destroy()
        
        self.SwapBuffers()

def doShapeGLPlot(mesh, angles, dirName):
    app = wx.PySimpleApp()
    frame = wx.Frame(None, wx.ID_ANY, "PCA GL Canvas", DEFAULT_POS, (400, 400))
    g = ShapeGLCanvas(frame, angles, dirName, mesh)
    frame.canvas = g
    frame.Show()
    app.MainLoop()
    app.Destroy()

def getUniformSphereAngles(level = 2):
    m = getSphereMesh(1, level)
    cosphi = m.VPos[:, 2]
    sinphi = np.sqrt(1-cosphi**2)
    costheta = m.VPos[:, 0]/sinphi
    sintheta = m.VPos[:, 1]/sinphi
    costheta[sinphi == 0] = 1
    sintheta[sinphi == 0] = 1
    cosphi[cosphi < -1] = -1
    cosphi[cosphi > 1] = 1
    costheta[costheta < -1] = -1
    costheta[costheta > 1] = 1
    phi = np.arccos(cosphi)
    theta = np.arccos(costheta)
    angles = np.zeros((len(phi), 2))
    angles[:, 0] = phi.flatten()
    angles[:, 1] = theta.flatten()
    angles = angles[np.argsort(angles[:, 0]), :]
    angles = angles[np.argsort(angles[:, 1]), :]
    plt.plot(angles[:, 0], 'b')
    plt.hold(True)
    plt.plot(angles[:, 1], 'r')
    plt.show()
    return angles

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Usage: python makeRotationImages.py <mesh name> <directory name> <NFrames>"
        sys.exit(0)
    [filename, dirName, NFrames] = [sys.argv[1], sys.argv[2], int(sys.argv[3])]
    np.random.seed(100) #For repeatable results randomly sampling
    
    if not os.path.exists(dirName):
        os.mkdir(dirName)
    
    m = PolyMesh()
    m.loadFile(filename)
    m.VPos = m.VPos*np.random.randn(m.VPos.shape[0], m.VPos.shape[1])
    
    #Output rotation video
    #angles = 2*np.pi*np.random.rand(NFrames, 2)
    
    #angles = np.zeros((NFrames, 2))
    #angles[:, 0] = np.linspace(0, 2*np.pi, NFrames+1)[0:NFrames]
    
    #angles = getUniformSphereAngles(3)
    
    [I, J] = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(0, 2*np.pi, 30))
    angles = np.zeros((I.size, 2))
    angles[:, 0] = I.flatten()
    angles[:, 1] = J.flatten()
    
    doShapeGLPlot(m, angles, dirName)


