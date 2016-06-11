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
    
    #Output rotation video
    angles = 2*np.pi*np.random.rand(NFrames, 2)
    doShapeGLPlot(m, angles, dirName)


