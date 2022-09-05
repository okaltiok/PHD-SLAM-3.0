function MonteCarloSimulations
    % the function can be used to generate data of N simulations, run the 
    % MCSs and compute the performance metrics averaged over the N MCSs
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 29/8/2022
    % Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.  

    % set option = 1 to create data
    % set option = 2 to run MCS
    % set option = 3 to compute performance metrics
    
    option = 1;
    
    % number of simulations
    N = 100;
    
    switch option
        case 1
            generate_data(N)
        case 2
            perform_mcs(N)
        case 3
            compute_performance_metrics
        otherwise
            fprintf('Incorrect option\n')
    end
end

function generate_data(N)
    reply = input('Are you certain, you want to generate new measurements? (y/n): ','s');
    if strcmpi(reply,'y')
        fprintf('generating new measurements\n')
    else
        fprintf('No action taken\n')
        return
    end

    % get current folder and define path names
    folderParts = strsplit(pwd, filesep);
    meas_path = strjoin([folderParts 'measurements\SD data set'], filesep);
    m_path = strjoin([folderParts(1:end-1) 'Shared Files\m source'], filesep);
    
    % add measurement path to measurements and m files
    if ~contains(path,meas_path), addpath(meas_path); end
    if ~contains(path,m_path), addpath(m_path); end
    
    % create data and save
    for ii = 1:N
        clear('params','model','sim','obj')
        [params,sim,filename]  = simulation_setup(ii,'save');
        save(filename, 'params', 'sim')
        fprintf('round %d\n',ii)
    end

    % remove created paths
    rmpath(meas_path);
    rmpath(m_path);
end

function perform_mcs(N)
    folderParts = strsplit(pwd, filesep);
    shared_m_path = strjoin([folderParts(1:end-1) 'Shared Files\shared m files'], filesep);
    if ~contains(path,shared_m_path), addpath(shared_m_path); end

    L = 5;
    J = 50;
    N_particle = logspace(0,2,3);
    resample = true;
    
    for n = N_particle
        POS_E = cell(N,1);
        THETA_E = cell(N,1);
        GOSPA_D = cell(N,1);
        CPU = cell(N,1);
        N_EFF = cell(N,1);

        pw = PoolWaitbar(N, 'Simulation in progress, please wait ...');
        
        % If Parallel Computing Toolbox isn't installed, replace
        % parfor-loop with a regular for-loop.
        parfor i = 1:N
            [pos_e,theta_e,gospa_d,cpu_time,Neff] = main(i,false,L,J,n,resample);

            POS_E{i,1} = pos_e;
            THETA_E{i,1} = theta_e;
            GOSPA_D{i,1} = gospa_d;
            CPU{i,1} = cpu_time;
            N_EFF{i,1} = Neff;

            increment(pw)
        end

        delete(pw)

        pos_e = cell2mat(POS_E);
        theta_e = cell2mat(THETA_E);
        gospa_d = cell2mat(GOSPA_D);
        cpu = cell2mat(CPU);

        fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], theta=%.2f [deg], gospa=%.2f [m], cpu=%.2f [ms]\n', ...
            J, ...
            L, ...
            n, ...
            sqrt(mean(pos_e.^2,'all','omitnan')), ...
            sqrt(mean(theta_e.^2,'all','omitnan'))*180/pi, ...
            mean(gospa_d(:,end),'all'), ...
            mean(cpu,'all')*1000);

        clearvars -except N_eff N POS_E THETA_E GOSPA_D CPU N_EFF J L N_particle resample i n

        filename = strcat('results/PHD_J',num2str(J),'L',num2str(L),'N',num2str(n),'.mat');
        save(filename,'POS_E','THETA_E','GOSPA_D','CPU','N_EFF')
    end
end

function compute_performance_metrics
    % define variables
    N_particle = logspace(0,2,3);
    L = 5;
    J = 50;
    
    % create filename
    filename = 'results/';
    file = [filename,'PHD_J%dL%dN%d.mat'];

    for n = N_particle
        % load data
        try
            load(sprintf(file,J,L,n));
        catch
            continue
        end
        pos_e = cell2mat(POS_E);
        theta_e = cell2mat(THETA_E);
        gospa_d = cell2mat(GOSPA_D);
        cpu = cell2mat(CPU);
        
        fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], theta=%.2f [deg], gospa=%.2f [m], cpu=%.2f [ms]\n', ...
            J, ...
            L, ...
            n, ...
            sqrt(mean(pos_e.^2,'all')), ...
            sqrt(mean(theta_e.^2,'all'))*180/pi, ...
            mean(gospa_d(:,end),'all'), ...
            mean(cpu,'all')*1000);
    end
end