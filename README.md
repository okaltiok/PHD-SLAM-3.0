1. Download the project to your computer.

2. PHD-SLAM 3.0 with synthetic data.\
    a) Open file ".\PHD-SLAM-3.0\Synthetic Data\MonteCarloSimulations.m"\
    b) Set option = 1 and run the code to create the data.\
    c) Set option = 2 and run the code to perform the Monte Carlo simulations (MCSs). This will take some time since it will perform 100 MCS with 1, 10 and 100        particles. You should get the following results:
    
       N:1, Pos.=3.45 [m], Head.=0.42 [deg], GOSPA=46.84 [m], ESS=100.00 [%]
       N:10, Pos.=2.66 [m], Head.=0.34 [deg], GOSPA=38.27 [m], ESS=49.34 [%]
       N:100, Pos.=2.18 [m], Head.=0.29 [deg], GOSPA=32.12 [m], ESS=49.21 [%]
    
    d) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3 and running the code. 
    
3. PHD-SLAM 3.0 with Victoria Park data set.\
    a) Open file ".\PHD-SLAM-3.0\Victoria Park\MonteCarloSimulations.m"\
    b) Set option = 1 and run the code to create the data.\
    c) Set option = 2 and run the code to perform MCSs. You should get the following result:
    
       N:1, pos=3.57 [m], ESS=100.00 [%]
    
    d) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3 and running the code. 
    
4. If Matlab throws an error from the mex-files, you need to compile the mex-files on your computer. The mex-files can be compiled using CompileCLibraries.m found in the folder ".\PHD-SLAM-3.0\Shared Files\mex source". Please note that the mex-files have to be compiled seperately for both data sets. Further instructions can be found in CompileCLibraries.m.     

5. The function perform_mcs() in MonteCarloSimulations.m utilizes Matlab's Parallel Computing Toolbox to run the simulations with multiple CPU cores in parallel. If you don't have Parallel Computing Toolbox installed, simply change the parfor-loop to a regular for-loop.
