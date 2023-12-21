1. Download the project to your computer.

2. PHD-SLAM 3.0 with synthetic data.\
    a) Open file ".\PHD-SLAM-3.0\Synthetic Data\MonteCarloSimulations.m"\
    b) Set option = 1 and run the code to create the data.\
    c) Set option = 2 and run the code to perform the Monte Carlo simulations (MCSs). This will take some time since it        will perform 100 MCS with 1, 10 and 100 particles. You should get the following results:
    
       N:1,   POS.=3.45 [m], HEAD.=0.42 [deg], GOSPA=46.84 [m], ESS=100.00 [%], RESAMP=0.00 [%]
       N:10,  POS.=2.66 [m], HEAD.=0.34 [deg], GOSPA=38.27 [m], ESS= 49.34 [%], RESAMP=5.14 [%]
       N:100, POS.=2.18 [m], HEAD.=0.29 [deg], GOSPA=32.12 [m], ESS= 49.21 [%], RESAMP=9.42 [%]
    
    d) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3         and running the code. 
    
3. PHD-SLAM 3.0 with Victoria Park data set.\
    a) Download "victoria_park.zip" from [1] and unzip the folder in to "Victoria Park/measurements/VP data set"
    b) Open file ".\PHD-SLAM-3.0\Victoria Park\MonteCarloSimulations.m"\
    c) Set option = 1 and run the code to create the data.\
    d) Set option = 2 and run the code to perform MCSs. You should get the following result:
    
       N:1,  POS.=3.57 [m], ESS=100.00 [%], RESAMP=0.00 [%]
       N:5,  POS.=3.38 [m], ESS= 25.17 [%], RESAMP=0.49 [%]
       N:10, POS.=3.36 [m], ESS= 13.88 [%], RESAMP=0.39 [%]
    
    e) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3         and running the code. 
    
5. If Matlab throws an error from the mex-files, you need to compile the mex-files on your computer. The mex-files can be compiled using CompileCLibraries.m found in the folder ".\PHD-SLAM-3.0\Shared Files\mex source". Please note that the mex-files have to be compiled seperately for both data sets. Further instructions can be found in CompileCLibraries.m.     

6. The function perform_mcs() in MonteCarloSimulations.m utilizes Matlab's Parallel Computing Toolbox to run the simulations with multiple CPU cores in parallel. If you don't have Parallel Computing Toolbox installed, simply change the parfor-loop to a regular for-loop.


References:

[1] “Victoria park SLAM data set,” Jose Guivant, Australian Centre for Field Robotics - The University of Sydney. Accessed Dec. 21, 2023. [Online]. Available: http://www-personal.acfr.usyd.edu.au/nebot/victoria_park.htm

