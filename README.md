# Project4

## project4 folder:

main.cpp takes four input arguments from the command line: 
- The first argument specifies the lattice size
- The second argument is the number of Monte Carlo cycles 
- The third argument is the temperature 
- The fourth argument is a variable that is set to either '1' or '2'. argv[4] == 1 tells the main program to loop over all Monte Carlo cycles. argv[4] == 2 tells the program to loop over the first 1000 Monte Carlo cycles and start sampling after an equilibriums state has been reached.

functions.cpp and functions.h contains functions that are called from main.cpp


## openmpi_project4 folder:

main.cpp runs with a custom executable, where we write '-n \<number of processors\> \<executable\>'. This 
specifies the number of processors we distribute the workload to for a given executable. 

functions_mpi.cpp and functions_mpi.h contains functions that are called from main.cpp


## files_python folder:

This folder contains output files and python codes used to obtain the results of the report. 
