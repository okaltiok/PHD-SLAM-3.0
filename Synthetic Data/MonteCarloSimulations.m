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
    
    option = 2;
    
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
        fprintf('generating new measurements\n')
    else
        fprintf('No action taken\n')
        return
    end

    folderParts = strsplit(pwd, filesep);
    meas_path = strjoin([folderParts 'measurements' 'SD data set'], filesep);
    m_path = strjoin([folderParts(1:end-1) 'Shared Files' 'm source'], filesep);
    
    % add measurement path to measurements and m files
    if ~contains(path,meas_path), addpath(meas_path); end
    if ~contains(path,m_path), addpath(m_path); end
    
    % create data and save
    for ii = 1:MCS
        clear('params','model','sim','obj')
        [params,sim,filename]  = simulation_setup(ii,'create');
        save(filename, 'params', 'sim')
        fprintf('round %d\n',ii)
    end

    % remove created paths
    rmpath(meas_path);
    rmpath(m_path);
end

function perform_mcs(MCS)
    folderParts = strsplit(pwd, filesep);
    shared_m_path = strjoin([folderParts(1:end-1) 'Shared Files' 'shared m files'], filesep);
    if ~contains(path,shared_m_path), addpath(shared_m_path); end

    % create filename
    if strcmp(filesep,'\')
         filename = strjoin([{'results\'} 'PHD_J%dL%dN%d.mat'], filesep);
    else
         filename = strjoin([{'results'} 'PHD_J%dL%dN%d.mat'], filesep);
    end
    
    % define variables
    N = 10;              % Particle number
    L = 5;              % Maximum number of IPL iterations
    J = 50;             % Maximum number of GM-OID components
    resample = true;    % resample flag
    
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
            main(i,false,L,J,N,resample);
        increment(pw)
    end
    
    delete(pw)
    
    pos_e = cell2mat(POS_E);
    theta_e = cell2mat(THETA_E);
    gospa_d = cell2mat(GOSPA_D);
    cpu = cell2mat(CPU);
    neff = cell2mat(N_EFF);
    
    fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], theta=%.2f [deg], gospa=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
        J, ...
        L, ...
        N, ...
        sqrt(mean(pos_e.^2,'all','omitnan')), ...
        sqrt(mean(theta_e.^2,'all','omitnan'))*180/pi, ...
        mean(gospa_d(:,end),'all'), ...
        mean(cpu,'all')*1000, ...
        mean(sum(cpu,2)), ...
        mean(neff,'all')*100);
    
    % save MCS results
    save(sprintf(filename,J,L,N),'POS_E','THETA_E','GOSPA_D','CPU','N_EFF')
end

function compute_performance_metrics
    % create filename
    if strcmp(filesep,'\')
        filename = strjoin([{'results\'} 'PHD_J%dL%dN%d.mat'], filesep);
    else
        filename = strjoin([{'results'}  'PHD_J%dL%dN%d.mat'], filesep);
    end
    
    % define variables
    N = 1;
    L = 5;
    J = 50;
    
    try
        load(sprintf(filename,J,L,N),'POS_E','THETA_E','GOSPA_D','CPU','N_EFF');
        
        pos_e = cell2mat(POS_E);
        theta_e = cell2mat(THETA_E);
        gospa_d = cell2mat(GOSPA_D);
        cpu = cell2mat(CPU);
        neff = cell2mat(N_EFF);
        
        fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], theta=%.2f [deg], gospa=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
            J, ...
            L, ...
            N, ...
            sqrt(mean(pos_e.^2,'all','omitnan')), ...
            sqrt(mean(theta_e.^2,'all','omitnan'))*180/pi, ...
            mean(gospa_d(:,end),'all'), ...
            mean(cpu,'all')*1000, ...
            mean(sum(cpu,2)), ...
            mean(neff,'all')*100);
    catch
        fprintf('Can''t load file, incorrect filename!\n')
    end
end