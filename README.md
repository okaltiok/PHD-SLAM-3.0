1. Download the project to your computer.

2. PHD-SLAM 3.0 with synthetic data.
    a) Open file ".\PHD-SLAM-3.0\Synthetic Data\MonteCarloSimulations.m"
    b) Set option = 1 and run the code to create MCS data. 
    c) Set option = 2 and run the code to perform MCS. This will take some time since it will perform 100 MCS using 1, 10 and 100 particles and the results            are stored in the results-folder.
    d) After the simulations are complete, set option = 3 and run the code to compute the performance metrics. You should get the following result:
    
J:50, L:5, N:1, pos=5.58 (4.25) [m], theta=0.63 (0.63) [deg], gospa=56.27 (40.61) [m], Neff=100.00 [%], Resampling=0.00 [%], cpu=0.35 (0.29) [ms]
J:50, L:5, N:10, pos=2.57 (1.70) [m], theta=0.32 (0.31) [deg], gospa=28.59 (17.87) [m], Neff=70.47 [%], Resampling=11.35 [%], cpu=2.51 (1.02) [ms]
J:50, L:5, N:100, pos=2.25 (1.48) [m], theta=0.29 (0.29) [deg], gospa=25.46 (15.69) [m], Neff=67.41 [%], Resampling=15.82 [%], cpu=24.38 (7.04) [ms]
    
3. PHD-SLAM 3.0 with Victoria Park data set.
    a) Open file ".\PHD-SLAM-3.0\Victoria Park\MonteCarloSimulations.m"
    b) Set option = 1 and run the code to create MCS data. 
    c) Set option = 2 and run the code to perform MCS. This will take some time since it will perform 100 MCS using 1, 10 and 100 particles and the results            are stored in the results-folder.
    d) After the simulations are complete, set option = 3 and run the code to compute the performance metrics. 
    
4. If Matlab throws an error from the mex-files, you need to possibly compile the mex-files. The mex-files can be compiled using CompileCLibraries.m and the mex-files have to be compiled seperately for both data sets. Further instructions can be found in CompileCLibraries.m.     

5. The function perform_mcs in MonteCarloSimulations.m utilizes Matlab's Parallel Computing Toolbox to run the simulations in parallel with multiple cores. If you don't have Parallel Computing Toolbox installed, simply change the parfor-loop to a regular for-loop.

6. 
