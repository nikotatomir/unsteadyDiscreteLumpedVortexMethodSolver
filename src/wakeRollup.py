import sys
import numpy as np

class wakeRollup():
    
    def __init__(self, wingSolution, wingKinematics, wakeConvection):
        
        self.wingSolution = wingSolution
        self.wingKinematics = wingKinematics
        self.wakeConvection = wakeConvection
        
        
        
        self.knownShedVortexCirculation = self.getKnownShedVortexCirculation()
        self.knownShedVortexPtEarth = self.getKnownShedVortexPtEarth()

    def getKnownShedVortexCirculation(self):
        
        knownShedVortexCirculation =  self.wingSolution.circulation[-1]
    
        return knownShedVortexCirculation

    def getKnownShedVortexPtEarth(self):
        
        if self.wakeConvection == "no":
                    
            knownShedVortexPtEarth = np.zeros(2)
            
            knownShedVortexPtEarth[:] = self.wingKinematics.unknownShedVortexPtEarth[:]
            
            return knownShedVortexPtEarth
        
        elif self.wakeConvection == "yes":
            print('\n')
            sys.exit('SYSTEM EXIT --> Error: Wake Convection Is Not Yet Implemented !!!')
        else:
            print('\n')
            sys.exit('SYSTEM EXIT --> Error: Wake Convection Must Be "yes" Or "no" !!!')
