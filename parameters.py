
import math
import numpy as np
import sys

#-------------------------------IMPORTANT--------------------------------------------#
# after the computations are done, there is an option to get an animation of the 
# airfoil motion. However, this takes a long time to compute (depending on the number)
# of time steps used. By default, it is turned of. If you want it, define it as 'yes' below
getAnimation = 'no' # 'yes' of 'no'

#------------------------BASIC-PROPERTIES--------------------------------------------#

c = 1 # chord length in meters
rho = 1 # density in kg/m^3

Mchord = 10 # number of panels
typeSpacing = "uniform" # "uniform" or "cosine"

# DO NOT CHANGE
wakeConvection = "no" # "yes" or "no", BUT "yes" option not implemented yet!

#------------------------MOTION DEFINITION-------------------------------------------#

# TRANSLATION
UinertialEarth = -50# m/s, x-component of intertial velocity defined in xe,ze system
# DO NOT CHANGE, WinertialEarth must be zero for panel Force to work
WinertialEarth = 0. # m/s,  z-component of intertial velocity defined in xe,ze system

# SINUSOIDAL PITCHING
pitchFreqHz = 12.732395*3/4 # in Hz
pitchAmpRad = math.radians(2) # value given in degrees, then converted to radians
pitchSSrad = math.radians(1) # value given in degrees, then converted to radians
rotationPt = 0.5*c # goes from 0 to c, DEFINED as "a" in Plots

#-------------------------TIME-PARAMETERS--------------------------------------------#

deltaT = 0.00075
NumberOfTimeSteps = 700
totalTime = NumberOfTimeSteps*deltaT 
time = np.arange(0, totalTime + 0.000001, deltaT)

#-------------------------FLOW-PARAMETERS--------------------------------------------#

Vmagnitude = np.linalg.norm([UinertialEarth,WinertialEarth])
pitchFreqRad = 2*np.pi*pitchFreqHz 
reducedFreq = (pitchFreqRad*c)/(2*Vmagnitude)
katzParameter = pitchAmpRad*pitchFreqRad*c/Vmagnitude
shedVortexFactor = 0.25 # should be between 0.2-0.3


#-------------------CONDITION-FOR-CONTINUATION-OF-SIMULATION-------------------------#
if (reducedFreq > 1) or (katzParameter > 1):
    print ("Reduced Frequency =", np.round(reducedFreq,3), "(-)", "// Should be Less Than 1")
    print ("Katz Parameter =", np.round(katzParameter,3), "(-)", "// Should be MUCH Less Than 1 \n")
    sys.exit("SYSTEM EXIT --> Error: reducedFreq is greater than 1 or katzParameter is greater than 1!!!")
else:
#-------------------PRINT-IMPORTANT-PARAMETERS-TO-SCREEN-----------------------------#
	print ("#----------------------------IMPORTANT FLOW PARAMETERS---------------------------------#")
	print ("Velocity Magnitude =", Vmagnitude,"(m/s) \n")
	print ("Steady State Angle of Attack =", np.round(math.degrees(pitchSSrad),3),"(in degrees)")
	print ("Amplitude of Angle of Attack Oscillations=", np.round(math.degrees(pitchAmpRad),3),"(in degrees)")
	print ("Pitch Frequency =", np.round(pitchFreqHz,3), "(in Hz)")
	if pitchFreqHz != 0:
		T = 1/pitchFreqHz # period for harmonical pitching
		numOscillations = totalTime/T
		print ("Pitch Period =", np.round(T,3), "(in seconds)")
		print ("Number of Oscillations Simulated =", np.round(numOscillations,3))
	else:
		print ("Pitch Period =", 0, "(in seconds)")
		print ("Number of Oscillations Simulated =", 0)
	print ("Reduced Frequency =", np.round(reducedFreq,3), "(-)", "// Should be Less Than 1")
	print ("Katz Parameter =", np.round(katzParameter,3), "(-)", "// Should be MUCH Less Than 1 \n")
	print ("Time Step =", np.round(deltaT,10), "(in seconds)")
	print ("Number of Time Steps =", NumberOfTimeSteps)
	print ("Final Time =", np.round(time[-1],10), "(in seconds)")
	print ("#--------------------------------------------------------------------------------------#\n")