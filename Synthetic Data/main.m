function varargout = main(varargin)
    % function main excecutes a single Monte Carlo simulation of PHD-SLAM 
    % using Synthetic data
    %
    % Input:
    %    varargin{1}         - file identifier
    %    varargin{2}         - plot mode
    %    varargin{3}         - N, particle number
    %    varargin{4}         - resample, true or false
    %
    % Output:
    %    varargout{1}        - (1 X T) vector of position errors
    %    varargout{2}        - (1 X T) vector of heading errors
    %    varargout{3}        - (1 X T) vector containing GOSPA
    %    varargout{4}        - (1 x T) vector of excecution times for each time instant
    %    varargout{5}        - (1 x T) vector that containts the N_eff for each time instant
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 29/9/2022
    % Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    % set to true to use MEX implementation
    MEX = false;
          
    % load simulation parameters
    if nargin == 0
        file_idx = 3;
        [params,sim,~] = simulation_setup(file_idx,'load',MEX);
        maxNumCompThreads(1); % force maximum number of computation threads to one
        clear plot_estimate
    else
        [params,sim,~] = simulation_setup(varargin{1},'load',MEX);
        params.f_mode = varargin{2}; 
        params.N_particle = varargin{3};
        params.resample = varargin{4};
        if params.resample
            params.N_eff = params.T_eff*params.N_particle;
        else
            params.N_eff = 0;
        end
    end

    % initialize PHD-SLAM density and arrays to store data
    [obj,sim] = initialize(sim,params);

    % last observation and control times
    t_prev = [];
    u_prev = [];
    
    % for each time instants, do
    for k = 1:params.T
        y = sim.y{k};
        
        % get control input and control times
        t = [t_prev sim.o_time{k}];
        u = [u_prev sim.odometry{k}];
        
        t_prev = sim.o_time{k};
        u_prev = sim.odometry{k};
        
        % perform one iteration of the PHD-SLAM filter
        [obj,est] = rbpf_phd(obj,y,u,t,params);

        % store data for evaluation purposes
        sim.XX(:,:,k) = [obj.xn];
        sim.WW(:,k) = [obj.w];
        sim.MM(:,k) = est.m_hat;
        sim.PP(:,:,k) = est.P_hat;
        sim.MM_map{1,k} = est.mu_hat;
        sim.PP_map{1,k} = est.C_hat;
        sim.Neff(1,k) = est.Neff;
        sim.cpu_time(k) = est.dt;
        
        % illustrate
        plot_estimate(obj,est,y,k,sim,params)
    end
    
    if nargin == 0
        params.f_mode = true;
    end
    [pos_e,theta_e,gospa_d] = performance_summary(params,sim);
    
    if nargin > 0 
        varargout{1} = pos_e;
        varargout{2} = theta_e;
        varargout{3} = gospa_d;
        varargout{4} = sim.cpu_time;
        varargout{5} = sim.Neff;
    end
    clear('params','sim','obj')
end