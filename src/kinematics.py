import numpy as np
import math
import sys
import matplotlib.pyplot as plt

class kinematics():
 
    def __init__(self, time, pitchAmpRad, pitchFreqHz, UinertialEarth, WinertialEarth, previousTEptEarth, shedVortexFactor, wingInitialGeometry):
        
        #-----------------------INPUT-PARAMETERS-------------------------------#
        
        self.previousTEptEarth = previousTEptEarth
        self.time = time
        self.pitchAmpRad = pitchAmpRad # (in radians)
        self.pitchFreqHz = pitchFreqHz # (in Hz)
        self.UinertialEarth = UinertialEarth 
        self.WinertialEarth = WinertialEarth
        self.shedVortexFactor = shedVortexFactor
        
        self.wingInitialGeometry = wingInitialGeometry
  
        #----------------------------------------------------------------------#
        
        self.pitchDisp = self.pitchAmpRad*math.sin(2*np.pi*self.pitchFreqHz*self.time)
        self.pitchVel = 2*np.pi*self.pitchFreqHz*self.pitchAmpRad*math.cos(2*np.pi*self.pitchFreqHz*self.time)
        self.translationDisp = np.array([self.UinertialEarth, self.WinertialEarth])*self.time
        self.translationVel = np.array([self.UinertialEarth, self.WinertialEarth])

        #----------------------------------------------------------------------#
        
        self.rotationMatrixEarthPt = np.array([[ math.cos(self.pitchDisp), math.sin(self.pitchDisp)],
                                             [-math.sin(self.pitchDisp), math.cos(self.pitchDisp)]])
        
        self.T_earthBody = np.array([[ math.cos(self.pitchDisp), -math.sin(self.pitchDisp)],
                                       [ math.sin(self.pitchDisp),  math.cos(self.pitchDisp)]])
    
        #----------------------------------------------------------------------#
        
        self.panelVortexPtsEarth = self.getPanelVortexPtsEarth()
        self.panelVortexPtsBody = self.getPanelVortexPtsBody()
        self.panelEvalPtsEarth = self.getPanelEvalPtsEarth()
        self.panelEvalPtsBody = self.getPanelEvalPtsBody()
        self.panelEndPtsEarth = self.getPanelEndPtsEarth()
                
        self.velocityPanelEvalPts_inertialEarth = self.getVelocityPanelEvalPts_inertialEarth()
        self.velocityPanelEvalPts_nonInertialBody = self.getVelocityPanelEvalPts_nonInertialBody()
        
        self.unknownShedVortexPtEarth = self.getUnknownShedVortexPtEarth()
        
    def getPanelVortexPtsEarth(self):
        panelVortexPtsEarth = np.zeros((self.wingInitialGeometry.Mchord,2))
        
        for i in range(self.wingInitialGeometry.Mchord):
            panelVortexPtsEarth[i,:] = self.translationDisp + np.matmul(self.rotationMatrixEarthPt, self.wingInitialGeometry.panelVortexPts[i,:]) 

        return panelVortexPtsEarth

    def getPanelVortexPtsBody(self):    
        panelVortexPtsBody = np.zeros((self.wingInitialGeometry.Mchord,2))
        
        originBodyInEarthAxes = self.translationDisp
        
        for i in range(self.wingInitialGeometry.Mchord):
            panelVortexPtsBody[i,:] = np.matmul(self.T_earthBody, self.panelVortexPtsEarth[i,:] - originBodyInEarthAxes)
        
        return panelVortexPtsBody
    
    def getPanelEvalPtsEarth(self):
        panelEvalPtsEarth = np.zeros((self.wingInitialGeometry.Mchord,2))
        
        for i in range(self.wingInitialGeometry.Mchord):
            panelEvalPtsEarth[i,:] = self.translationDisp + np.matmul(self.rotationMatrixEarthPt, self.wingInitialGeometry.panelEvalPts[i,:]) 

        return panelEvalPtsEarth

    def getPanelEvalPtsBody(self):    
        panelEvalPtsBody = np.zeros((self.wingInitialGeometry.Mchord,2))
        
        originBodyInEarthAxes = self.translationDisp
        
        for i in range(self.wingInitialGeometry.Mchord):
            panelEvalPtsBody[i,:] = np.matmul(self.T_earthBody, self.panelEvalPtsEarth[i,:] - originBodyInEarthAxes)
        
        return panelEvalPtsBody

    def getPanelEndPtsEarth(self):
        panelEndPtsEarth = np.zeros((self.wingInitialGeometry.Mchord+1,2))
        
        for i in range(self.wingInitialGeometry.Mchord+1):
            panelEndPtsEarth[i,:] = self.translationDisp + np.matmul(self.rotationMatrixEarthPt, self.wingInitialGeometry.panelEndPts[i,:]) 

        return panelEndPtsEarth

    def getVelocityPanelEvalPts_inertialEarth(self):
        velocityPanelEvalPts_inertialEarth = np.zeros((self.wingInitialGeometry.Mchord, 2))
        
        originBodyInEarthAxes = self.translationDisp
        
        for i in range(self.wingInitialGeometry.Mchord):
            rVecEarth = self.panelEvalPtsEarth[i,:]-originBodyInEarthAxes
            velocityPanelEvalPts_inertialEarth[i,:] = ( self.translationVel + np.array([self.pitchVel*rVecEarth[1], -self.pitchVel*rVecEarth[0] ]) ) 
        
        return velocityPanelEvalPts_inertialEarth

    def getVelocityPanelEvalPts_nonInertialBody(self):
        velocityPanelEvalPts_nonInertialBody = np.zeros((self.wingInitialGeometry.Mchord, 2))
        
        for i in range(self.wingInitialGeometry.Mchord):
            velocityPanelEvalPts_nonInertialBody[i,:] = np.matmul(self.T_earthBody, -self.velocityPanelEvalPts_inertialEarth[i,:])
                
        return velocityPanelEvalPts_nonInertialBody

    def getUnknownShedVortexPtEarth(self):
        unknownShedVortexPtEarth = np.zeros(2)
        currentTEptEarth = self.panelEndPtsEarth[-1,:]

        unknownShedVortexPtEarth[:] = ( self.shedVortexFactor*(self.previousTEptEarth - currentTEptEarth) ) + currentTEptEarth
        
        return unknownShedVortexPtEarth


