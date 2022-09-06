1. Download the project to your computer.

2. PHD-SLAM 3.0 with synthetic data.\
    a) Open file ".\PHD-SLAM-3.0\Synthetic Data\MonteCarloSimulations.m"\
    b) Set option = 1 and run the code to create the data.\
    c) Set option = 2 and run the code to perform the Monte Carlo simulations (MCSs). This will take some time since it will perform 100 MCS with 1, 10 and 100 particles. You should get the following result:
    
            J:50, L:5, N:1, pos=5.58 [m], theta=0.63 [deg], gospa=56.27 [m], cpu=0.27 [ms]
            J:50, L:5, N:10, pos=2.57 [m], theta=0.32 [deg], gospa=28.59 [m], cpu=2.29 [ms]
            J:50, L:5, N:100, pos=2.25 [m], theta=0.29 [deg], gospa=25.46 [m], cpu=23.10 [ms]
    
    d) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3 and running the code. 
    
3. PHD-SLAM 3.0 with Victoria Park data set.\
    a) Open file ".\PHD-SLAM-3.0\Victoria Park\MonteCarloSimulations.m"\
    b) Set option = 1 and run the code to create the data.\
    c) Set option = 2 and run the code to perform MCSs. This will take some time since it will perform 100 MCS with 1, 10 and 100 particles. You should get the following result:
    
            J:50, L:5, N:1, pos=4.37 [m], cpu=0.58 [ms]
            J:50, L:5, N:10, pos=2.84 [m], cpu=4.49 [ms]
            J:50, L:5, N:100, pos=2.88 [m], cpu=43.17 [ms]
    
    d) The results are stored in the results-folder, and the performance metrics can be computed by setting option = 3 and running the code. 
    
4. If Matlab throws an error from the mex-files, you need to possibly compile the mex-files. The mex-files can be compiled using CompileCLibraries.m and the mex-files have to be compiled seperately for both data sets. Further instructions can be found in CompileCLibraries.m.     

5. The function perform_mcs() in MonteCarloSimulations.m utilizes Matlab's Parallel Computing Toolbox to run the simulations with multiple CPU cores in parallel. If you don't have Parallel Computing Toolbox installed, simply change the parfor-loop to a regular for-loop.
