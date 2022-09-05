function varargout = main(varargin)
    % function main excecutes a single Monte Carlo simulation of PHD-SLAM using
    % the Victoria Park data set
    %
    % Input:
    %    varargin{1}         - file identifier
    %    varargin{2}         - plot mode
    %    varargin{3}         - L, max number of OID iterations
    %    varargin{4}         - J, max number of GMM OID components
    %    varargin{5}         - N, particle number
    %    
    % Output:
    %    varargout{1}        - (1 X T) vector that containts the RMSE of vehicle for each time instant
    %    varargout{2}        - (1 x T) vector that containts the CPU time for each time instant
    %    varargout{3}        - (1 x T) vector that containts the N_eff for each time instant
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 24/8/2022
    % Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    
    % load simulation parameters
    if nargin == 0
        [params,sim,~] = simulation_setup(1,'load');
        maxNumCompThreads(1); % force maximum number of computation threads to one
        clear plot_estimate
    else
        [params,sim,~] = simulation_setup(varargin{1},'load');
        params.f_mode = varargin{2}; 
        params.L = varargin{3};
        params.J = varargin{4};
        params.N_particle = varargin{5};
    end 
     
    % get current folder and define path names
    folderParts = strsplit(pwd, filesep);
    mex_path = strjoin([folderParts 'compiled_mex'], filesep);
    m_path = strjoin([folderParts(1:end-1) 'Shared Files\m source'], filesep);
    shared_m_path = strjoin([folderParts(1:end-1) 'Shared Files\shared m files'], filesep);
    
    % add / remove paths
    if ~contains(path,shared_m_path), addpath(shared_m_path); end
    if params.MEX
        if ~contains(path,mex_path), addpath(mex_path); end
        if contains(path,m_path), rmpath(m_path); end
    else
        if ~contains(path,m_path), addpath(m_path); end
        if contains(path,mex_path), rmpath(mex_path); end
    end
    
    % initialize the random number generator
    rng(params.seed)

    params.T = size(sim.u_time,2);
    params.K = size(sim.o_time,2);
    
    % initialize arrays to store data
    sim.XX = zeros(params.xn_dim,params.N_particle,params.K);
    sim.WW = zeros(params.N_particle,params.K);
    sim.MM = zeros(params.xn_dim,params.K);
    sim.PP = zeros(params.xn_dim,params.xn_dim,params.K);
    sim.MM_map = cell(1,params.K);
    sim.PP_map = cell(1,params.K);
    sim.cpu_time = zeros(1,params.K);
    sim.Neff = zeros(1,params.K); 
    
    % initialize PHD-SLAM density
    obj = initialize(sim,params);

    % help variables
    k_last_control = 1;
    j = 1;
    t_prev = [];
    u_prev = [];
     
    % for each time instants, do
    for k = 1:params.T   
        if sim.o_time{j} < sim.u_time{k}

            y = sim.y{j};
           
            % get control input and control times
            t = [t_prev sim.u_time{k_last_control:k-1} sim.o_time{j}];
            u = [u_prev sim.odometry{k_last_control:k-1}];
            if k > 1
                t_prev = sim.o_time{j};
                u_prev = sim.odometry{k-1};
            end
            
            % perform one iteration of the filter
            [obj,est] = rbpf_phd(obj,y,u,t,params);
             
            % store data for evaluation purposes
            sim.XX(:,:,j) = [obj.xn];
            sim.WW(:,j) = [obj.w];
            sim.MM(:,j) = est.m_hat;
            sim.PP(:,:,j) = est.P_hat;
            sim.MM_map{1,j} = est.mu_hat;
            sim.PP_map{1,j} = est.C_hat;
            sim.Neff(1,j) = est.Neff;
            sim.cpu_time(j) = est.dt;
            
            % illustrate 
            plot_estimate(obj,est,y,j,sim,params)
            
            k_last_control = k;
            j = j + 1;
        end

        if j > params.K 
            break
        end
    end
    
    [sim,pos_e] = performance_summary(params,sim);

    if nargin > 0 
        varargout{1} = pos_e;
        varargout{2} = sim.cpu_time;
        varargout{3} = sim.Neff;
    end
    clear('params','sim','obj','map')
end