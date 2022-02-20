import numpy as np
import math
import sys
import matplotlib.pyplot as plt

class initialGeometry():
 
    def __init__(self, c, Mchord, typeSpacing, pitchSSrad, rotationPt):
        self.c = c
        self.Mchord = Mchord
        self.typeSpacing = typeSpacing
        self.pitchSSrad = pitchSSrad
        self.rotationPt = rotationPt
                
        self.panelEndPts = self.getPanelEndPts()    
        self.panelVortexPts = self.getPanelVortexPts()
        self.panelEvalPts = self.getPanelEvalPts()
        self.panelNormVec = self.getPanelNormVec()
        self.panelTangVec = self.getPanelTangVec()
        self.panelArea = self.getPanelArea()
        self.checkPanelNormVectorOrthogonality = self.getCheckPanelNormVectorOrthogonality()
        self.checkPanelNormVectorMagnitude = self.getCheckPanelNormVectorMagnitude()
               
    def getPanelEndPts(self):
        panelEndPts = np.zeros((self.Mchord+1,2))
        thetaSpacing = np.linspace(0,180,self.Mchord+1)
        rotationMatrix = np.array([[math.cos(self.pitchSSrad), math.sin(self.pitchSSrad)],[-math.sin(self.pitchSSrad), math.cos(self.pitchSSrad)]])
        
        # Computation of panelEndPts x-component
        if self.typeSpacing == "uniform":
            panelEndPts[:,0] = np.linspace(0,self.c,self.Mchord+1) - self.rotationPt
        elif self.typeSpacing == "cosine":
            for i in range(self.Mchord+1):
               panelEndPts[i,0] = 0.5*self.c*(1-math.cos(math.radians(thetaSpacing[i]))) - self.rotationPt
        else:
            sys.exit('typeSpacing must be "uniform" or "cosine"')

        # Rotation of panelEndPts by pitchSS angle
        for i in range(self.Mchord+1):
            panelEndPts[i,:] = np.matmul(rotationMatrix, panelEndPts[i,:])
            
        return panelEndPts
    
    def getPanelVortexPts(self):
        panelVortexPts = np.zeros((self.Mchord,2))
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:]
            panelVortexPts[i,:] = self.panelEndPts[i,:] + 0.25*vec
        
        return panelVortexPts
    
    def getPanelEvalPts(self):
        panelEvalPts = np.zeros((self.Mchord,2))
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:]
            panelEvalPts[i,:] = self.panelEndPts[i,:] + 0.75*vec
        
        return panelEvalPts
    
    def getPanelNormVec(self):
        panelNormVec = np.zeros((self.Mchord,2))
        
        thetaRotation = np.pi/2
        transformationMatrix = np.array([[math.cos(thetaRotation), -math.sin(thetaRotation)],[math.sin(thetaRotation), math.cos(thetaRotation)]])
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:] # this vector is rotated by pi/2 in the counter-clockwise direction
            panelNormVec[i,:] = ( np.matmul(transformationMatrix, vec) ) / np.linalg.norm(vec)
        return panelNormVec

    def getPanelTangVec(self):
        panelTangVec = np.zeros((self.Mchord,2))
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:]
            panelTangVec[i,:] = vec / np.linalg.norm(vec)
        
        return panelTangVec

    def getPanelArea(self):
        panelArea = np.zeros(self.Mchord)
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:]
            panelArea[i] = np.linalg.norm(vec)
            
        return panelArea
    
    def getCheckPanelNormVectorOrthogonality(self):
        checkPanelNormVectorOrthogonality = np.zeros(self.Mchord)
        
        for i in range(self.Mchord):
            vec = self.panelEndPts[i+1,:] - self.panelEndPts[i,:] 
            checkPanelNormVectorOrthogonality[i] = np.vdot(vec,self.panelNormVec[i,:]) 
            
        return checkPanelNormVectorOrthogonality
        
    def getCheckPanelNormVectorMagnitude(self):
        checkPanelNormVectorMagnitude = np.zeros(self.Mchord)
        
        for i in range(self.Mchord):
            checkPanelNormVectorMagnitude[i] = np.linalg.norm(self.panelNormVec[i,:])
        
        return checkPanelNormVectorMagnitude

    def gridlines(self):
        plt.minorticks_on()
        plt.grid(zorder = 0, which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75')
  
    def getInitialGeometryPlot(self):
        fig = plt.figure(0)
        ax = fig.add_subplot(111, aspect='equal', autoscale_on = False, xlim=(-1.5*self.c,1.5*self.c), ylim=(-1.5*self.c,1.5*self.c))
        ax.plot(self.panelEndPts[:,0], self.panelEndPts[:,1], 'r-', linewidth = 0.5, label = 'Thin Wing (airfoil)')
        ax.plot(self.panelEndPts[:,0], self.panelEndPts[:,1], 'r|', markersize = 5, label = 'Panel End Points' )
        ax.plot(self.panelVortexPts[:,0], self.panelVortexPts[:,1], 'b.', markersize = 2, label = 'Bound Vortex Points')
        ax.plot(self.panelEvalPts[:,0], self.panelEvalPts[:,1], 'c.', markersize = 2, label = 'Evaluation Points')
        ax.quiver(self.panelEvalPts[:,0], self.panelEvalPts[:,1], self.panelNormVec[:,0], self.panelNormVec[:,1], scale = 20, width = 0.003)
        ax.legend()
        self.gridlines()
        plt.savefig("airfoilGeometry_t=0", bbox_inches='tight', dpi = 250)
        print ("\nInitial Geometry At t=0s Plotted")