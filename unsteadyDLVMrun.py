from src.initialGeometry import initialGeometry
from src.kinematics import kinematics
from src.solution import solution
from src.systemProperties import systemProperties
from src.wakeRollup import wakeRollup

import postProcess

#---------------------------SIMULATION-PARAMETERS------------------------------# 
from parameters import *

#-------------OBJECTS-FOR-MONITORING-SIMULATION-RESULTS------------------------# 
wingKinematics = np.empty(len(time), dtype = 'object')
wingSystemProperties = np.empty(len(time), dtype = 'object')
wingWakeRollup = np.empty(len(time), dtype = 'object')

#-------------------INITIAL-AIRFOIL-POSITION,time t=0--------------------------#
#note: at this time, the non-inertial body coordinate system xb,zb and the inertial earth coordinate system xe,ze overlap
        
wingInitialGeometry = initialGeometry(c, Mchord, typeSpacing, pitchSSrad, rotationPt)
wingInitialGeometry.getInitialGeometryPlot()
        
TElocation = wingInitialGeometry.panelEndPts[-1,:]

print ("\nSimulation Started")
print ("\nTime Step: ", end="")
for index in range(1,len(time)):
    print (index, end=", ", flush=True)

    #---------------------------KINEMATICS-----------------------------------------#
    wingKinematics[index] = kinematics(time[index], pitchAmpRad, pitchFreqHz, UinertialEarth, WinertialEarth, TElocation, shedVortexFactor, wingInitialGeometry)
    # updating TE location
    TElocation = wingKinematics[index].panelEndPtsEarth[-1,:]
            
    #---------------------------SOLUTION--------------------------------------------#
    wingSolution = solution(wingInitialGeometry, wingKinematics[index], index, deltaT, wingSystemProperties[index-1], wingWakeRollup[0:index])
    
    #------------------------SYSTEM-PROPERTIES--------------------------------------#
    wingSystemProperties[index] = systemProperties(wingSolution, wingKinematics[index], wingInitialGeometry, wingWakeRollup[0:index], index, wingSystemProperties[index-1], rho, Vmagnitude, deltaT)
    
    #---------------------------WAKE-ROLLUP-----------------------------------------#
    # adding position of latest shed free vortex
    wingWakeRollup[index] = wakeRollup(wingSolution, wingKinematics[index], wakeConvection)
    
print ("\n\nSimulation Done")


#---------------------------POST-PROCESSING------------------------------# 

if int(pitchFreqHz) == 0:
    postProcess.impulsiveFlatPlatePlot(wingSystemProperties)
else:
    postProcess.harmonicOscillationPlot(wingKinematics, wingSystemProperties)

if getAnimation == 'yes':
    postProcess.animationMovieKinematics(wingKinematics)
elif getAnimation == 'no':
    pass
else:
    print('\n')
    sys.exit('SYSTEM EXIT --> Error: getAnimation Must Be "yes" Or "no" !!!')


