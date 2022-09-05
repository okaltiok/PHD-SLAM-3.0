function [params,sim,filename] = simulation_setup(file_idx,mode)
    % This function is used to generate the data and configure the 
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
    % Last Rev : 1/9/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    % if the data doesn't exist in the workspace
    if ~(evalin('base','exist(''params'',''var'')') && ...
            evalin('base','exist(''sim'',''var'')') && ...
            evalin('base','exist(''filename'',''var'')'))

        switch mode
            % load data
            case 'load'
                j = 0;
                while 1
                    try
                        filename = ['.\measurements\measurement_' num2str(file_idx)];
                        load(filename), break;
                    catch
                        j = j + 1;
                        fprintf('load error: %d\n',j);
                        pause(1)
                        if j >= 10
                            return
                        end
                    end
                end
                
            % generate data and save
            case 'save'
                % create filename
                filename = ['.\measurements\measurement_' num2str(file_idx)];
                
                % initialize parameters
                params = initialize_parameters;
                
                % initialize random seed based on the file identifier
                rng(file_idx,'twister');
                seed = rng;
                params.seed = seed;
                
                % load map, robot trajectory and robot controls
                load('Environment_SLAM.mat', ...
                    'Landmark_Groundtruth', ...
                    'Robot_Groundtruth', ...
                    'Robot_Control');
                
                map = Landmark_Groundtruth(:,2:3)'; % landmarks
                
                sim = '';
                
                % initial state and arrays
                sim.state = zeros(params.xn_dim,params.T);
                sim.P0 = zeros(params.xn_dim);
                sim.x0 = Robot_Groundtruth(:,1);
                sim.odometry = cell(1,params.T);
                sim.u_time = cell(1,params.T);
                sim.y = cell(1,params.T);
                sim.o_time = cell(1,params.T);
                sim.map = map;
                
                % generate trajectory and measurements
                x = sim.x0;
                t = 0;
                sqrtQ = chol(params.Qu,'lower');
                sqrtR = chol(params.R,'lower');
                for k = 1:params.T
                    if k > 1
                        % propagate vehicle
                        x = f_func(x,Robot_Control(:,k-1),params.dt);
                    end
                    
                    % Select set of landmarks that are visible within vehicle's semi-circular field-of-view
                    dx = map(1,:) - x(1);
                    dy = map(2,:) - x(2);
                    range2 = dx.^2 + dy.^2;
                    theta = mod(atan2(dy,dx)-x(3) + pi,2*pi) - pi;
                    
                    idx = range2 <= params.fov_range.^2 & abs(theta) <= params.fov_angle;
                    
                    lm = map(:,idx);
                    
                    % number of landmarks and clutter measurements
                    N = size(lm,2);
                    N_clutter = poissrnd(params.lambda,1);
                    
                    % create noiseless measurements
                    h = zeros(params.h_dim,N+N_clutter);
                    for n = 1:N+N_clutter
                        if n <= N
                            h(:,n) = h_func(x,lm(:,n));
                        else
                            random_range = rand*params.fov_range;
                            random_angle = rand*(params.fov_angle*2) - params.fov_angle;
                            random_location = [x(1) + random_range*cos(random_angle + x(3)); ...
                                               x(2) + random_range*sin(random_angle + x(3))];
                            h(:,n) = h_func(x,random_location);
                        end
                    end
                    
                    % add noise and remove meas. that are misdetected
                    y = h + sqrtR*randn(params.h_dim,N+N_clutter);
                    y = y(:,rand(1,N+N_clutter) <= params.P_D);
                    
                    % handle wrap around
                    y(2,:) = mod(y(2,:) + pi,2*pi) - pi;

                    % store measurements 
                    sim.state(:,k) = x;
                    sim.y{1,k} = y;
                    sim.o_time{1,k} = t;
                    sim.odometry{1,k} = Robot_Control(:,k) + sqrtQ*randn(params.u_dim,1);
                    sim.u_time{1,k} = t;
                    
                    t = t + params.dt;
                end
        end
    else
        params = evalin('base','params');
        sim = evalin('base','sim');
        filename = evalin('base','filename');
    end
end


function params = initialize_parameters
    % This method initializes the parameters

    % general
    params.MEX = true;                                             % MEX flag
    params.f_mode = true;                                          % plot flag
    params.visualize_passive_components = true;                    % plot components outside FOV flag
    params.T = 4000;                                               % simulation length
    params.dt = 1;                                                 % sampling interval
    
    % FOV parameters
    params.fov_range = 150;                                        % FoV range
    params.fov_angle = 90*pi/180;                                  % FoV angle
    params.range_buffer = 10;                                      % parameter used to determine if a landmark near the FOV
    params.angle_buffer = 10*pi/180;                               % parameter used to determine if a landmark near the FOV

    % PHD parameters
    params.lambda = 5;                                             % expectation of the Poisson distribution
    params.lambda_c = params.lambda/ ...
        ((2*params.fov_angle/(2*pi))*pi*params.fov_range^2);       % clutter intensity
    params.P_D  = 0.95;                                            % probability of detection
    params.P_B  = 1e-6;                                            % birth probability
    params.P_G = 1e-9;                                             % gating tail probability
    params.gating_size = chi2inv(1 - params.P_G,2);                % gating threshold value
    params.w_min =  log(1e-6);                                     % pruning threshold
    params.merging_threshold = 50;                                 % merging threshold
    params.eta_threshold = 0.7;                                    % threshold of landmark estimate

    % PF and OID parameters
    params.resample = false;                                       % resample flag
    params.N_particle = 1;                                         % number of particles
    params.N_eff = 0;                                              % effective sample size
    params.J = 50;                                                 % max number of assignments computed using Murty's algorithm
    params.DA_threshold = log(10^(-3));                            % threshold to stop Murty's algorithm at n-th iteration, if gain[1] - gain[n] < params.DA_threshold 
    params.L = 5;                                                  % max number of OID iterations
    params.epsilon = 1e-3;                                         % threshold value used to evaluate convergence of the iterative OID approximation
    params.kappa = 1e-3;                                           % test statistic tail probability  
    params.gamma = chi2inv(1 - params.kappa , (0:100)*2);          % test statistic value used to evaluate goodness of the linearization
    
    % model parameters
    params.xn_dim = 3;                                             % vehicle dimension
    params.xl_dim = 2;                                             % landmark dimension                 
    params.u_dim = 2;                                              % control dimension
    params.h_dim = 2;                                              % measurement dimension
    params.Qn = diag([1e-9 1e-9 1e-18]);                           % Covariance of vehicle state (x, y, heading)
    params.Qu = diag([0.8^2 (0.5*pi/180)^2]);                      % Covariance of control input (velocity, steering)
    params.R = diag([0.8^2 (0.3*pi/180)^2]);                       % Covariance of measurements  (range, angle)
end



