1. Download project to your computer.

2. PHD-SLAM 3.0 with synthetic data.
    a) Open file ".\PHD-SLAM-3.0\Synthetic Data\MonteCarloSimulations.m"
    b) Set option = 1 and run the code to create MCS data. 
    c) Set option = 2 and run the code to perform MCS. This will take some time since it will perform 100 MCS using 1, 10 and 100 particles and the results            are stored in the results-folder.
    d) After the simulations are complete, set option = 3 and run the code to compute the performance metrics. 
    
3. PHD-SLAM 3.0 with Victoria Park data set.
    a) Open file ".\PHD-SLAM-3.0\Victoria Park\MonteCarloSimulations.m"
    b) Set option = 1 and run the code to create MCS data. 
    c) Set option = 2 and run the code to perform MCS. This will take some time since it will perform 100 MCS using 1, 10 and 100 particles and the results            are stored in the results-folder.
    d) After the simulations are complete, set option = 3 and run the code to compute the performance metrics. 
    
4. If Matlab throws an error from the mex-files you need to check the compiler settings and possibly compile the mex-files.     
