README

NOTE: boost and ODEINT are not included in this distribution and need to be downloaded seperately as explained below

Files:

compile_mex_files.m:
Script to help compile the mex files in this project, see note about Boost for required dependencies

time_comparison.m
Script to compare the time it takes to compute solutions to the system of ODEs using the three implementations in this project (mexode, mexodestiff, matlabs ode15s)

mexodestiff1.cpp:
A mex function that uses the Rosenbrock method to solve the system of ODEs

mexodestiff1.hpp:
Header file for above.


COMPILING THE MEX FILES
The mex files all rely on boost and ODEINT

1. Download Boost:
http://www.boost.org/users/download/

2. Download ODEINT:
http://headmyshoulder.github.com/odeint-v2/downloads.html

3. Extract Boost, and extract ODEINT

4. Copy the folder "boost" from ODEINT's extracted directory to the folder in BOOST that already contains a folder called boost

5. Allow the folders to be merged

6. Use the provided compilation script from Matlab

7. To compile by hand from matlab use:
mex -I"/path/to/boost" mexodestiff1.cpp