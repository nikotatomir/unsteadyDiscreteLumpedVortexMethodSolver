import math
import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from parameters import *

def gridlines():
    plt.minorticks_on()
    plt.grid(zorder = 0, which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
    plt.grid(zorder = 0, which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75')
    plt.grid(zorder = 0, which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    plt.grid(zorder = 0, which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75')

def impulsiveFlatPlatePlot(wingSystemProperties):
    print ("\nPlotting Unsteady Lift & Circulation Of Impulsively Started Flat Plate...")

    nonDimTimeList = []
    nonDimTotalPanelCirculationList = []
    nonDimTotalPanelLiftCoefficientList = []
    
    wagnerList = []
    
    totalLiftCoefficient_SS = 2*np.pi*pitchSSrad
    totalCirculation_SS = np.pi*pitchSSrad*c*abs(UinertialEarth)

    for i in range(1,len(time)):

        nonDimTimeList.append(2*time[i]*abs(UinertialEarth)/c)
        wagnerList.append( 1 - (0.165*math.exp(-0.045*nonDimTimeList[i-1])) - (0.335*math.exp(-0.3*nonDimTimeList[i-1])) )
        
        nonDimTotalPanelCirculationList.append( wingSystemProperties[i].totalPanelCirculation / totalCirculation_SS )
        nonDimTotalPanelLiftCoefficientList.append( wingSystemProperties[i].totalPanelLiftCoefficient / totalLiftCoefficient_SS )
        
    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    ax.plot(nonDimTimeList, wagnerList, 'k-', linewidth = 0.5, label = '$Wagner$ $Function$')
    ax.plot(nonDimTimeList, nonDimTotalPanelCirculationList, 'b-', linewidth = 0.5, label = '$\\frac{\Gamma}{\Gamma_{\infty}}$')
    ax.plot(nonDimTimeList, nonDimTotalPanelLiftCoefficientList, 'm-', linewidth = 0.5, label = '$\\frac{c_l}{c_{l,\infty}}$')
    gridlines()
    ax.set_xlabel("Non-Dimensional Time (2Ut/c)")
    ax.set_ylabel("Non-Dimensional Circulation $\\frac{\Gamma}{\Gamma_{\infty}}$ & Lift Coefficient $\\frac{c_l}{c_{l,\infty}}$")
    ax.set_title("Impulsively Started Airfoil")
    ax.legend()
    ax.set_xlim(0, nonDimTimeList[-1])
    ax.set_ylim(0,1.2)
    
    plt.savefig("unsteadyLiftAndCirculation.png", bbox_inches='tight', dpi = 250)
    print ("\nUnsteady Lift & Circulation Of Impulsively Started Flat Plate Plotted")

def harmonicOscillationPlot(wingKinematics, wingSystemProperties):
    print ("\nPlotting Unsteady Lift Of Harmonically Oscillating Airfoil...")

    theodorsenFunc = scipy.special.hankel2(1,reducedFreq) / ( scipy.special.hankel2(1,reducedFreq) + ( 1j*scipy.special.hankel2(0,reducedFreq) ) )
    a = ( rotationPt - (c/2) ) / (c/2)
    clalpha_SS = 2*np.pi
    theoCoeff = ( 1j*np.pi*reducedFreq ) + (a*np.pi*(reducedFreq**2)) + (clalpha_SS*theodorsenFunc) + ( (clalpha_SS*theodorsenFunc)*(1j*reducedFreq*(0.5-a)) )
    
    
    aoaList = []
    totalPanelLiftCoefficientList = []
    
    theodorsenLiftReal = []
    theodorsenAOAreal= []
    
    totalLiftCoefficient_SS = []
    
    period = 1/pitchFreqHz
    lastPeriodIndexStart = -period/deltaT - 2
    
    for i in range(int(lastPeriodIndexStart),0):
        
        aoaList.append( math.degrees(wingKinematics[i].pitchDisp + pitchSSrad))
        totalPanelLiftCoefficientList.append( wingSystemProperties[i].totalPanelLiftCoefficient ) 
        
        theodorsenLift = (clalpha_SS*pitchSSrad) + (theoCoeff*pitchAmpRad*np.exp(1j*pitchFreqRad*time[i]))
        theodorsenAOA = pitchSSrad + ( pitchAmpRad*np.exp(1j*pitchFreqRad*time[i]) )
        
        theodorsenLiftReal.append( np.real(theodorsenLift) )
        theodorsenAOAreal.append( math.degrees( np.real(theodorsenAOA) ) )
    
        totalLiftCoefficient_SS.append( 2*np.pi*(wingKinematics[i].pitchDisp + pitchSSrad)) 
    
    fig = plt.figure(4)
    ax = fig.add_subplot(111)
    ax.plot(aoaList, totalPanelLiftCoefficientList, 'm-', linewidth = 0.5, label = '$Unsteady$ $Lift$ $Coefficient$')
    ax.plot(theodorsenAOAreal, theodorsenLiftReal, 'r-', linewidth = 0.5, label = '$Unsteady$ $Theodorsen$ $Lift$ $Coefficient$')
    ax.plot(aoaList, totalLiftCoefficient_SS, 'b-', linewidth = 0.5, label = '$Steady$ $Lift$ $Coefficient$' )
    gridlines()
    ax.set_xlabel("Angle of Attack $\\theta$ (in degrees)")
    ax.set_ylabel("Lift Coefficient $c_l$ (-)")
    ax.set_title(f"Harmonically Oscillating Airfoil\nReduced Freq. $k={np.round(reducedFreq,3)}$, Rotation Pt $a={rotationPt}c$, Chord $c={c}m$")
    ax.legend()
    ax.set_xlim(np.round(min(aoaList)), np.round(max(aoaList)))
    #ax.set_ylim(0.0,0.25)
    
    plt.savefig("unsteadyLift.png", bbox_inches='tight', dpi = 250)
    print ("\nUnsteady Lift Of Harmonically Oscillating Airfoil Plotted")

def animationMovieKinematics(wingKinematics):
    print ("\nCreating Animation...")
    # set up the figure and subplot

    fig = plt.figure(2,figsize=(32.0,12.0))
    fig.canvas.set_window_title('Matplotlib Animation')
    ax = fig.add_subplot(111, aspect='equal', autoscale_on = False, xlim=(-15,2), ylim=(-2,2))
    ax.grid()
    ax.set_title('2D Wing Motion Animation', fontsize = 22)
    ax.set_xlabel('Distance (m)', fontsize = 18)
    line, = ax.plot([],[], 'r|-', lw=1, label = 'Panel End Points') 
    line2, = ax.plot([],[], 'b.', markersize=8, label = 'Bound Vortex Points') 
    line3, = ax.plot([],[], 'c.', markersize=8, label = 'Evaluation Points') 
    line4, = ax.plot([],[], 'g.', markersize=2, label = 'Free Vortex Points') 
    ax.legend(fontsize = 18)

    
    def init():
        line.set_data([],[])
        line2.set_data([],[])
        line3.set_data([],[])
        line4.set_data([],[])
        
        return line, line2, line3, line4, 
    
    def animate(i):
        #print (i)

        x_points = wingKinematics[i+1].panelEndPtsEarth[:,0]
        z_points =  wingKinematics[i+1].panelEndPtsEarth[:,1]
        line.set_data(x_points, z_points)
        
        x_points = wingKinematics[i+1].panelVortexPtsEarth[:,0]
        z_points = wingKinematics[i+1].panelVortexPtsEarth[:,1]                 
        line2.set_data(x_points, z_points)
        
        x_points = wingKinematics[i+1].panelEvalPtsEarth[:,0]
        z_points = wingKinematics[i+1].panelEvalPtsEarth[:,1]                 
        line3.set_data(x_points, z_points)
        
        x_points = []
        z_points = []
        for j in range(i+1):
            x_points.append(wingKinematics[j+1].unknownShedVortexPtEarth[0])
            z_points.append(wingKinematics[j+1].unknownShedVortexPtEarth[1])
        line4.set_data(x_points, z_points)
        
        return line, line2, line3, line4,
    
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(wingKinematics)-1, interval=1, blit=True, repeat=False)
    ani.save('airfoilMotionAnimation.gif', fps = 30)
    print ("\nAnimation Created")