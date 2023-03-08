function [params,sim,filename] = simulation_setup(file_idx,mode,MEX)
    % This function is used to generate the data and configures the 
    % parameters or load the data from a file
    %
    % Input:
    %    file_idx   - file identifier
    %    mode       - either 'load' or 'save'
    %
    % Output:
    %    params     - simulation parameters
    %    sim        - struct containing the simulation data
    %    filename   - file identifier
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 3/10/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    % create filename
    filename = strjoin([{'.'}, 'measurements', ...
        strcat('measurement_',num2str(file_idx))], filesep);
    
    switch mode
        case 'load'
            [params, sim] = load_data(filename,MEX);
        case 'create'
            [params, sim] = create_data(file_idx);
        otherwise
            fprintf('Incorrect option!\n');
    end
end

function [params, sim] = load_data(filename,MEX)
    % load data
    j = 0;
    while 1
        try
            loaded_struct = load(filename,'params','sim');
            params = loaded_struct.params;
            sim = loaded_struct.sim;
            break;
        catch
            j = j + 1;
            fprintf('load error: %d\n',j);
            pause(1)
            if j >= 10
                return
            end
        end
    end
       
    % store MEX flag
    params.MEX = MEX;
    
    % define Matlab-paths, get current folder and define path names
    folderParts = strsplit(pwd, filesep);
    shared_m_path = strjoin([folderParts(1:end-1) ['Shared Files',filesep,'shared m files']], filesep);
    mex_path = strjoin([folderParts 'compiled_mex'], filesep);
    m_path = strjoin([folderParts(1:end-1) ['Shared Files',filesep,'m source']], filesep);
    
    % add / remove paths
    if ~contains(path,shared_m_path), addpath(shared_m_path); end
    if params.MEX
        if ~contains(path,mex_path), addpath(mex_path); end
        if contains(path,m_path), rmpath(m_path); end
    else
        if ~contains(path,m_path), addpath(m_path); end
        if contains(path,mex_path), rmpath(mex_path); end
    end
    
    % remove path to other data set
    synthetic_data_path = strjoin([folderParts(1:end-1) ['Synthetic Data',filesep,'compiled_mex']], filesep);
    if contains(path,synthetic_data_path), rmpath(synthetic_data_path); end

    % initialize the random number generator
    rng(params.seed)    
end

function [params, sim] = create_data(file_idx)
    % initialize parameters
    params = initialize_parameters;

    % initialize random seed based on the file identifier
    rng(file_idx,'twister');
    seed = rng;
    params.seed = seed;

    % load ground truth
    load('aa3_gpsx.mat','La_m','Lo_m','timeGps');
    sim.state = [Lo_m La_m]';
    sim.gps_time = timeGps'/1000;

    % load trajectory and controls
    load('aa3_dr.mat','speed','steering','time')

    % initial state
    sim.x0 = [-67.6493; -41.7142; 35.5*pi/180];
    sim.P0 = zeros(3);

    % store control information
    T = size(time,1);
    sim.odometry = cell(1,T);
    sim.u_time = cell(1,T);
    for k = 1:T
        u = [speed(k,1) steering(k,1)]';
        sim.odometry{1,k} = u;
        sim.u_time{1,k} = time(k,1)/1000;
    end

    % load measurements
    load('aa3_lsr2.mat','LASER','TLsr')

    K = size(TLsr,1);
    sim.y = cell(1,K);
    sim.o_time = cell(1,K);

    global AAr;
    AAr = (0:360)*pi/360 ;

    % store measurements and append random clutter to meas.
    j = 1;
    for k = 1:K
        % get measurements from the tree detector
        observations = detectTreesI16(double(LASER(k,:))/100);

        if ~isempty(observations)
            % convert from 0-pi to -pi/2 - pi/2
            observations(2,:) = observations(2,:) - pi/2;

            % observations within FoV
            idx = observations(1,:) <=  params.fov_range & ...
                abs(observations(2,:)) <=  params.fov_angle & ...
                observations(3,:) >= 0.0;

            y_detections = observations(1:2,idx);

            % create clutter
            if params.clutter
                N_clutter = poissrnd(params.lambda,1);
                random_range = rand(1,N_clutter)*params.fov_range;
                random_angle = rand(1,N_clutter)*(params.fov_angle*2) - params.fov_angle;
                y_clutter = [random_range; random_angle];
            else
                y_clutter = [];
            end

            % append measurements to cell array
            sim.y{1,j} = [y_detections y_clutter];
            sim.o_time{1,j} = double(TLsr(k,1))/1000;
            j = j + 1;
        end
    end
    sim.o_time(j:end) = [];
    sim.y(j:end) = [];
    
    params.T = size(sim.o_time,2);
    params.K = size(sim.u_time,2);
end
function params = initialize_parameters
    % This method initializes the parameters

    % general
    params.MEX = true;                                             % MEX flag
    params.f_mode = true;                                          % plot flag
    params.visualize_passive_components = true;                    % plot components outside FOV flag
    params.clutter = true;                                         % clutter measurements flag
    
    % FOV parameters
    params.fov_range = 50;                                         % FoV range
    params.fov_angle = 85*pi/180;                                  % FoV angle
    params.range_buffer = 5;                                       % parameter used to determine if a landmark near the FOV 
    params.angle_buffer = 5*pi/180;                                % parameter used to determine if a landmark near the FOV 

    % PHD parameters
    params.lambda = 5;                                             % expectation of the Poisson distribution
    params.lambda_c = params.lambda/ ...
        ((2*params.fov_angle/(2*pi))*pi*params.fov_range^2);       % clutter intensity
    params.P_D  = 0.7;                                             % probability of detection
    params.P_B  = 1e-6;                                            % birth probability
    params.P_G = 1e-9;                                             % gating tail probability 
    params.gating_size = chi2inv(1 - params.P_G,2);                % gating threshold value
    params.w_min =  log(1e-6);                                     % pruning threshold
    params.merging_threshold = 10;                                 % merging threshold
    params.etaT = log((1-params.P_D)^2);                           % threshold of landmark estimate

    % PF and OID parameters
    params.resample = false;                                       % resample flag
    params.N_particle = 1;                                         % number of particles
    params.T_eff = 0.2;                                            % resampling threshold
    params.N_eff = 0;                                              % effective sample size
    params.J = 50;                                                % max number of assignments computed using Murty's algorithm
    params.DA_threshold = log(10^(-3));                            % threshold to stop Murty's algorithm at n-th iteration, if gain[1] - gain[n] < params.DA_threshold 
    params.L = 5;                                                  % max number of OID iterations
    params.epsilon = 1e-3;                                         % threshold value used to evaluate convergence of the iterative OID approximation

    % model parameters
    params.xn_dim = 3;                                             % vehicle dimension
    params.xl_dim = 2;                                             % landmark dimension                 
    params.u_dim = 2;                                              % control dimension
    params.h_dim = 2;                                              % measurement dimension
    params.Qn = diag([1e-9 1e-9 1e-18]);                           % Covariance of vehicle state (x, y, heading)
    params.Qu = diag([1^2 (4*pi/180)^2]);                          % Covariance of control input (velocity, steering)
    params.R = diag([1^2 (1*pi/180)^2]);                           % Covariance of measurements  (range, angle)    
end



