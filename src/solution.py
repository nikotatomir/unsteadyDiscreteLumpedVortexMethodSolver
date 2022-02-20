import numpy as np

class solution():
    
    def __init__(self, wingInitialGeometry, wingKinematics, index, deltaT, previousWingSystemProperties, wingWakeRollup):
        
        self.wingInitialGeometry = wingInitialGeometry
        self.wingKinematics = wingKinematics  
        self.index = index
        self.deltaT = deltaT
        self.previousWingSystemProperties = previousWingSystemProperties
        self.wingWakeRollup = wingWakeRollup
        
        self.influenceCoefficientMatrix = self.getInfluenceCoefficientMatrix()
        self.RHS = self.getRHS()
        self.circulation = self.getCirculation()

    def unitStrengthVortexInducedVelocity(self, generalEvalPt, generalVortexPt):
        
        vec = generalEvalPt- generalVortexPt
        r = np.linalg.norm(vec)
                
        unitStrengthVortexInducedVelocity = ( 1/(2*np.pi*r*r) ) * np.matmul(np.array([[0,1],[-1,0]]), vec)
        
        return unitStrengthVortexInducedVelocity
    
    def getInfluenceCoefficientMatrix(self):
        influenceCoefficientMatrix = np.zeros( ( self.wingInitialGeometry.Mchord+1 , self.wingInitialGeometry.Mchord+1 ) )
        
        # COMPUTING INFLUENCE COEFFICIENT OF BOUND VORTICES
        for i in range(self.wingInitialGeometry.Mchord):
            for j in range(self.wingInitialGeometry.Mchord):
                
                unitStrengthPanelVortexInducedVelocity_nonInertialBody = self.unitStrengthVortexInducedVelocity( self.wingInitialGeometry.panelEvalPts[i,:], 
                                                                                                                 self.wingInitialGeometry.panelVortexPts[j,:] ) 
                
                influenceCoefficientMatrix[i,j] = np.vdot( unitStrengthPanelVortexInducedVelocity_nonInertialBody , 
                                                           self.wingInitialGeometry.panelNormVec[i,:] )
        
        # COMPUTING INFLUENCE COEFFICIENT OF ONE UNKNOWN SHED VORTEX
        for i in range(self.wingInitialGeometry.Mchord):
            
            unitStrengthShedVortexInducedVelocity_inertialEarth =  self.unitStrengthVortexInducedVelocity( self.wingKinematics.panelEvalPtsEarth[i,:], 
                                                                                                           self.wingKinematics.unknownShedVortexPtEarth[:] ) 
            
            unitStrengthShedVortexInducedVelocity_nonInertialBody = np.matmul( self.wingKinematics.T_earthBody, 
                                                                                unitStrengthShedVortexInducedVelocity_inertialEarth )
            
            influenceCoefficientMatrix[i,-1] =  np.vdot( unitStrengthShedVortexInducedVelocity_nonInertialBody , 
                                                         self.wingInitialGeometry.panelNormVec[i,:] )
        
        influenceCoefficientMatrix[-1,:] = 1        
        
        return influenceCoefficientMatrix

    def getRHS(self):
        RHS = np.zeros(self.wingInitialGeometry.Mchord+1)
        
        if self.index == 1:

            for i in range(self.wingInitialGeometry.Mchord):
                RHS[i] = -np.vdot(self.wingKinematics.velocityPanelEvalPts_nonInertialBody[i,:], self.wingInitialGeometry.panelNormVec[i,:])
            
            # NOTE that RHS[-1] = 0 automatically because RHS = np.zeros(  )
            
        else:

            knownShedVortexInducedVelocity_nonInertialBody = np.zeros( ( self.wingInitialGeometry.Mchord, self.index-1 , 2 ) )
            
            for i in range(self.wingInitialGeometry.Mchord):
                for j in range(1, self.index): # loops over all the known shed vortices from previous time steps
                    knownShedVortexInducedVelocity_inertialBody = ( self.wingWakeRollup[j].knownShedVortexCirculation*self.unitStrengthVortexInducedVelocity( self.wingKinematics.panelEvalPtsEarth[i,:], 
                                                                                                                                                          self.wingWakeRollup[j].knownShedVortexPtEarth[:]) ) 
                  
                    knownShedVortexInducedVelocity_nonInertialBody[i,j-1,:] =  np.matmul(self.wingKinematics.T_earthBody, knownShedVortexInducedVelocity_inertialBody)
            

                RHS[i] = - ( np.vdot(self.wingKinematics.velocityPanelEvalPts_nonInertialBody[i,:] , self.wingInitialGeometry.panelNormVec[i,:]) +
                             np.vdot( (sum(knownShedVortexInducedVelocity_nonInertialBody[i]))     , self.wingInitialGeometry.panelNormVec[i,:] ) )
            
            RHS[-1] = self.previousWingSystemProperties.totalPanelCirculation
            
        return RHS

    def getCirculation(self):
        circulation = np.zeros(self.wingInitialGeometry.Mchord+1)
        
        circulation =  np.matmul(np.linalg.inv(self.influenceCoefficientMatrix), self.RHS)
        
        return circulation