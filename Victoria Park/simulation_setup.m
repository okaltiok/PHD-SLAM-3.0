function [params,sim,filename] = simulation_setup(file_idx,mode)
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
                        filename = ['.\measurements\journal_measurements\measurement_' num2str(file_idx)];
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
                
            % generate data and sace
            case 'save'
                % create filename
                filename = ['.\measurements\measurement_' num2str(file_idx)];
                
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
    params.eta_threshold = 0.7;                                    % threshold of landmark estimate

    % PF and OID parameters
    params.resample = false;                                       % resample flag
    params.N_particle = 1;                                         % number of particles
    params.N_eff = 0;                                              % effective sample size
    params.J = 50;                                                % max number of assignments computed using Murty's algorithm
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
    params.Qu = diag([1^2 (4*pi/180)^2]);                          % Covariance of control input (velocity, steering)
    params.R = diag([1^2 (1*pi/180)^2]);                           % Covariance of measurements  (range, angle)    
end



