function MonteCarloSimulations
    % the function can be used to generate data of MC simulations, run the 
    % MCSs and compute the performance metrics averaged over the MCSs
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

    % set option = 1 to create data
    % set option = 2 to run MCS
    % set option = 3 to compute performance metrics
    
    option = 1;
    
    % number of simulations
    MCS = 100;
    
    switch option
        case 1
            generate_data(MCS)
        case 2
            perform_mcs(MCS)
        case 3
            compute_performance_metrics
        otherwise
            fprintf('Incorrect option\n')
    end
end

function generate_data(MCS)
    reply = input('Are you certain, you want to generate new measurements? (y/n): ','s');
    if strcmpi(reply,'y')
        fprintf('Generating data for Monte Carlo simulations.\n')
    else
        fprintf('No action taken.\n')
        return
    end

    folderParts = strsplit(pwd, filesep);
    meas_path = strjoin([folderParts 'measurements' 'SD data set'], filesep);
    m_path = strjoin([folderParts(1:end-1) 'Shared Files' 'm source'], filesep);
    shared_path = strjoin([folderParts(1:end-1) 'Shared Files' 'shared m files'], filesep);
    
    % add measurement path to measurements and m files
    if ~contains(path,meas_path), addpath(meas_path); end
    if ~contains(path,m_path), addpath(m_path); end
    if ~contains(path,shared_path), addpath(shared_path); end
    
    % create data and save
    pw = PoolWaitbar(MCS, 'Creating data, please wait ...');
    for ii = 1:MCS
        clear('params','model','sim','obj')
        [params,sim,filename]  = simulation_setup(ii,'create');
        save(filename, 'params', 'sim')
        increment(pw)
    end
    delete(pw)

    % remove created paths
    rmpath(meas_path);
    rmpath(m_path);
    rmpath(shared_path);
end

function perform_mcs(MCS)
    folderParts = strsplit(pwd, filesep);
    shared_m_path = strjoin([folderParts(1:end-1) 'Shared Files' 'shared m files'], filesep);
    if ~contains(path,shared_m_path), addpath(shared_m_path); end

    % create filename
    if strcmp(filesep,'\')
         filename = strjoin([{'results\'} 'PHD_N%d.mat'], filesep);
    else
         filename = strjoin([{'results'} 'PHD_N%d.mat'], filesep);
    end
    
    N = logspace(0,2,3)';    % Particle number
    resample = true;        % resample flag
    for j = 1:size(N,1)
        POS_E = cell(MCS,1);
        THETA_E = cell(MCS,1);
        GOSPA_D = cell(MCS,1);
        CPU = cell(MCS,1);
        N_EFF = cell(MCS,1);

        % If Parallel Computing Toolbox isn't installed, replace
        % parfor-loop with a regular for-loop.
        pw = PoolWaitbar(MCS, 'Simulation in progress, please wait ...');
        parfor i = 1:MCS
            [POS_E{i,1},THETA_E{i,1},GOSPA_D{i,1},CPU{i,1},N_EFF{i,1}] = ...
                main(i,false,N(j),resample);
            increment(pw)
        end
        delete(pw)
        
        gospa_d = cell2mat(GOSPA_D);
        fprintf('N:%d, Pos.=%.2f [m], Head.=%.2f [deg], GOSPA=%.2f [m], ESS=%.2f [%%]\n', ...
            N(j), ...
            sqrt(mean(cell2mat(POS_E).^2,'all','omitnan')), ...
            sqrt(mean(cell2mat(THETA_E).^2,'all','omitnan'))*180/pi, ...
            mean(gospa_d(:,end),'all'), ...
            mean(cell2mat(N_EFF),'all')*100);
        
        % save MCS results
        save(sprintf(filename,N(j)),'POS_E','THETA_E','GOSPA_D','CPU','N_EFF')
    end
end

function compute_performance_metrics
    % create filename
    if strcmp(filesep,'\')
        filename = strjoin([{'results\'} 'PHD_N%d.mat'], filesep);
    else
        filename = strjoin([{'results'} 'PHD_N%d.mat'], filesep);
    end
    
    % Particle number
    N = logspace(0,2,3)';    
    for j = 1:size(N,1)
        try
            load(sprintf(filename,N(j)),'POS_E','THETA_E','GOSPA_D','CPU','N_EFF');
            
            gospa_d = cell2mat(GOSPA_D);
            fprintf('N:%d, Pos.=%.2f [m], Head.=%.2f [deg], GOSPA=%.2f [m], ESS=%.2f [%%]\n', ...
                N(j), ...
                sqrt(mean(cell2mat(POS_E).^2,'all','omitnan')), ...
                sqrt(mean(cell2mat(THETA_E).^2,'all','omitnan'))*180/pi, ...
                mean(gospa_d(:,end),'all'), ...
                mean(cell2mat(N_EFF),'all')*100);
        catch
            fprintf('Can''t load file, incorrect filename!\n')
        end
    end
end
