The plate size is 0.02x0.3*0.5m with the following resolutions.
plate2.obj v 578 f 1152
plate2a.obj v 2306 f 4608

The code is described as follows, 
result.txt...the value of Kirchhoff tensor
plate1.obj...plate2.obj
plate2.obj...plate2a.obj
readobj.m....transform .obj file to .mat format
TestPlate.m...execute computation of Kirchhoff tensor
Kirchhoff3D.m...computate Kirchhoff tensor
simulator.m...the Lie-group rigid body integrator
DrawTraj.m...rendering the simulation results


I tried to fix out the computation error of the approach in contrast to the following work 
Equation (D3) [1] pp.26
Equation (6.4) [2] pp.16
The comparison error is not antcipated.


[1]W. B. Wang et.al., Influence of aspect ratio on tumbling plates,J. Fluid Mech. (2013)
[2]A. Anderson et.al., Unsteady aerodynamics of fluttering and tumbling plates, J. Fluid Mech. (2005)

--
HAORAN XIE