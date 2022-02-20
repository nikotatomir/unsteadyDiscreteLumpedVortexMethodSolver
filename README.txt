This code uses the unsteady discrete lumped-vortex panel method (DLVM) algorithm described in Katz & Plotkin in Section 13.10, pg. 407 [1] to compute the aerodynamic response to the following unsteady thin-airfoil motions:

	1. harmonically pitching motion about the steady state angle of attack.
	2. impulsively started thin airfoil motion at a steady state angle of attack.

The harmonically pitching motion is validated with Theodorsen Theory.
The impulsively started thin airfoil motion is validated with the Wagner Function. 
Some examples of the code results along with the settings of the paramaters.py file for the specific case are given in the sampleResults/ folder. 

The simulation solutions are outputed in 2 plots & 1 video:
	1. A plot showing the airfoil geometry at time is 0s.
	2. A plot showing the unsteady lift.
	3. A video showing the airfoil motion.

The code is run by executing the "unsteadyDLVMrun.py" file. Nothing in this file needs to be altered. The parameters for the simulation are defined in the "parameters.py" file.

A a few remarks about the "parameters.py" file:
	Here the getAnimation variable is by default 'no'. This is because if you increase the number of airfoil panels (i.e. finer discretization), the rendering time for the video increases substansially. My advice is, just for the video, to decrease the number of panels to 2 by setting the Mchord=2. The airfoil motion will still be the same, it will just take less time to render. 
	The pitchFreqHz variable defines whether the motion is harmoincally pitching or an impusively started thin airfoil. If it is nonzero it is the former, if it is zero it is the latter.
	Another important thing to note is that when running the simulation for a harmonically pitching airfoil, the "totalTime" variable should be larger than the period of 1 oscillation - otherwise the "postProcess.py" script will not work.  
	Finally, you can run the "parameters.py" file separately and it will print to the screen the input parameters and the parameters calculated from the inputs. In this way, you can double check your simulations conditions before running the main code.

NOTE: 
	Only thin symmetrical airfoils are considered; however, the code can be easily extended to include thin cambered airfoils.
	The wake rollup is not fully implemented, i.e. the effect of the induced velocities of the system vortices on the wake is not taken into account. However, the results provided in the sampleResults/ folder show that the the aerodynamic loads still correlate well with the theory.
	The blockScheme.pdf file gives an idea of how the instances and classes in the code interact with each other.
	
Package requirements for running the simulations in a Conda virtual environment --> requirements.yml file.

Bibliography

[1] Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics. Vol. 13. Cambridge university press, 2001.

Code written in 2021.