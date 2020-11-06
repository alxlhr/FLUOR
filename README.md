# FLUOR
Underwater ray tracing code developed in Python during my internship

How to run :

Compile the fortran subroutines with :
f2py -c -m fortran_loss fortran_loss.f90 
/
f2py -c -m fortran_arr fortran_arr.f90

Edit the user.py file to your needs, and run it.


Capacities : 
- 1D and 2D sound speed fields
- Variable bathymetry
- Transmission Loss
- Arrival calculations

Problems :
- Bad user interface
- Overestimation of the direct ray
- Discrepancies with Bellhop (is it really a problem ?)
- No makefile
