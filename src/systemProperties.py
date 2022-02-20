import numpy as np
import math as math

class systemProperties():
    
    def __init__(self, wingSolution, wingKinematics, wingInitialGeometry, wingWakeRollup, index, previousWingSystemProperties, rho, Vmagnitude, deltaT):
        
        self.wingSolution = wingSolution
        self.wingKinematics = wingKinematics
        self.wingInitialGeometry = wingInitialGeometry
        self.wingWakeRollup = wingWakeRollup
        self.index = index
        self.rho = rho
        self.Vmagnitude = Vmagnitude
        self.deltaT = deltaT
        
        if self.index > 1:
            self.previousWingSystemProperties = previousWingSystemProperties
        
        # CHECK THIS MATRIX
        self.T_bodyEarth = np.array([[  math.cos(self.wingKinematics.pitchDisp),  math.sin(self.wingKinematics.pitchDisp)],
                                     [ -math.sin(self.wingKinematics.pitchDisp),  math.cos(self.wingKinematics.pitchDisp)]])
        
        self.panelCirculation = self.getPanelCirculation()
        self.totalPanelCirculation = self.getTotalPanelCirculation()
        
        self.wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody = self.getWakeInducedVelocitiesAtPanelEvalPts_nonInertialBody()
        self.deltaPotential = self.getDeltaPotential()
        self.panelDeltaPressure = self.getPanelDeltaPressure()
        self.panelLiftForceEarth = self.getPanelLiftForceEarth()
        self.totalPanelLiftForceEarth = self.getTotalPanelLiftForceEarth()
        self.totalPanelLiftCoefficient = self.getTotalPanelLiftCoeffient()

        
    def getPanelCirculation(self):
        panelCirculation = np.zeros(self.wingInitialGeometry.Mchord)
        
        panelCirculation[:] = self.wingSolution.circulation[0:-1]
        
        return panelCirculation

    def getTotalPanelCirculation(self):
        
        totalPanelCirculation = sum ( self.wingSolution.circulation[0:-1] )
        
        return totalPanelCirculation

# PRESSURE COMPUTATIONS

    def getWakeInducedVelocitiesAtPanelEvalPts_nonInertialBody(self):
        
        wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody = np.zeros( ( self.wingInitialGeometry.Mchord, self.index, 2) )
        
        latestKnownShedVortexCirculation =  self.wingSolution.circulation[-1]
        latestKnownShedVortexPtEarth =  self.wingKinematics.unknownShedVortexPtEarth[:]
        
        
        if self.index == 1:
            # velocity induced by latest shed vortex on panelEvalPts
            for i in range( self.wingInitialGeometry.Mchord ):
                wakeInducedVelocitiesAtPanelEvalPts_inertialBody = latestKnownShedVortexCirculation*self.wingSolution.unitStrengthVortexInducedVelocity( self.wingKinematics.panelEvalPtsEarth[i,:], 
                                                                                                                                                         latestKnownShedVortexPtEarth[:]) 
                
                wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody[i,self.index-1,:] = np.matmul(self.wingKinematics.T_earthBody, wakeInducedVelocitiesAtPanelEvalPts_inertialBody)
        else:
            # velocity induced ALL known shed vortices from previous time steps on panelEvalPts
            for i in range( self.wingInitialGeometry.Mchord ):
                for j in range (1,self.index): # number of known shed vortices ( = ALL known shed vortices from previous time steps)
                    #print("Known Shed Vortex Number=", j)
                    wakeInducedVelocitiesAtPanelEvalPts_inertialBody = ( self.wingWakeRollup[j].knownShedVortexCirculation*self.wingSolution.unitStrengthVortexInducedVelocity( self.wingKinematics.panelEvalPtsEarth[i,:], 
                                                                                                                                                                                self.wingWakeRollup[j].knownShedVortexPtEarth[:]) ) 
                      
                    wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody[i,j-1,:] =   np.matmul(self.wingKinematics.T_earthBody, wakeInducedVelocitiesAtPanelEvalPts_inertialBody)
                    
            # velocity induced by latest shed vortex on panelEvalPts
            for i in range( self.wingInitialGeometry.Mchord ):
                wakeInducedVelocitiesAtPanelEvalPts_inertialBody = latestKnownShedVortexCirculation*self.wingSolution.unitStrengthVortexInducedVelocity( self.wingKinematics.panelEvalPtsEarth[i,:], 
                                                                                                                                                         latestKnownShedVortexPtEarth[:]) 
                
                wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody[i,self.index-1,:] = np.matmul(self.wingKinematics.T_earthBody, wakeInducedVelocitiesAtPanelEvalPts_inertialBody)
        
        return wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody

    def getDeltaPotential(self):
        
        deltaPotential = np.zeros(self.wingInitialGeometry.Mchord)
        
        for i in range( self.wingInitialGeometry.Mchord ):
            deltaPotential[i] = sum(self.panelCirculation[0:i+1])

        return deltaPotential    
    
    def getPanelDeltaPressure(self):

        # for index == 1, wingWakeRollup is empty, for index >1, use wingWakeRollup[1:index]
        
        panelDeltaPressure = np.zeros(self.wingInitialGeometry.Mchord)
        
        if self.index == 1:
            
            for i in range(self.wingInitialGeometry.Mchord):
                term1 = np.vdot( self.wingKinematics.velocityPanelEvalPts_nonInertialBody[i,:] + sum(self.wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody[i]) , self.wingInitialGeometry.panelTangVec[i,:] )
                term2 = self.panelCirculation[i] / self.wingInitialGeometry.panelArea[i]
                term3 = self.deltaPotential[i] / self.deltaT
                
                panelDeltaPressure[i] = self.rho*( (term1*term2) + term3 )
            
        else:
            for i in range(self.wingInitialGeometry.Mchord):
                term1 = np.vdot( self.wingKinematics.velocityPanelEvalPts_nonInertialBody[i,:] + sum(self.wakeInducedVelocitiesAtPanelEvalPts_nonInertialBody[i]) , self.wingInitialGeometry.panelTangVec[i,:] )
                term2 = self.panelCirculation[i] / self.wingInitialGeometry.panelArea[i]
                term3 = ( self.deltaPotential[i] - self.previousWingSystemProperties.deltaPotential[i] ) / self.deltaT
                
                panelDeltaPressure[i] = self.rho*( (term1*term2) + term3 )
                
        return panelDeltaPressure
        

    def getPanelLiftForceEarth(self):
        
        #print ("panel Force is incorrect if WinertialEarth is nonzero!")
        
        panelLiftForceEarth = np.zeros((self.wingInitialGeometry.Mchord,2))
        
        for i in range( self.wingInitialGeometry.Mchord ):
            #should there be a minus sign before self.panelDeltaPressure ...
            panelLiftForceBody = ( self.panelDeltaPressure[i]*self.wingInitialGeometry.panelArea[i] )*self.wingInitialGeometry.panelNormVec[i,:]
            panelLiftForceEarth[i,:] = np.matmul(self.T_bodyEarth, panelLiftForceBody)
        
        return panelLiftForceEarth

    def getPanelLiftCoefficient(self):
        pass

    def getTotalPanelLiftForceEarth(self):
        
        totalPanelLiftForceEarth = sum(self.panelLiftForceEarth)
        
        return totalPanelLiftForceEarth

    def getTotalPanelLiftCoeffient(self):
        
        totalPanelLiftCoefficient = self.totalPanelLiftForceEarth[1] / (0.5*self.rho*self.Vmagnitude*self.Vmagnitude*self.wingInitialGeometry.c)
        
        return totalPanelLiftCoefficient

    def getPanelMoment(self): # about which point?
        pass
    
    
    

    