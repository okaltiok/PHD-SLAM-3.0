function MonteCarloSimulations
    % the function can be used to generate data of N simulations, run the 
    % MCSs and compute the performance metrics averaged over the N MCSs
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 3/10/2022
    % Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.  

    % set option = 1 to create data
    % set option = 2 to run MCS
    % set option = 3 to compute performance metrics
    
    option = 3;
    
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

    % get current folder, define path name and add it to the path
    folderParts = strsplit(pwd, filesep);
    meas_path = strjoin([folderParts 'measurements' 'VP data set'], filesep);
    if ~contains(path,meas_path), addpath(meas_path); end

    % create data and save
    for ii = 1:MCS
        clear('params','model','sim','obj')
        [params,sim,filename]  = simulation_setup(ii,'create');
        save(filename, 'params', 'sim')
        fprintf('round %d\n',ii)
    end
    
    % remove created paths
    rmpath(meas_path);
end

function perform_mcs(MCS)
    % get current folder, define path name and add it to the path
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
    N = 1;  % Particle number
    L = 5;  % Maximum number of IPL iterations
    J = 50; % Maximum number of GM-OID components

    POS_E = cell(MCS,1);
    CPU = cell(MCS,1);
    N_EFF = cell(MCS,1);
    
    % If Parallel Computing Toolbox isn't installed, replace
    % parfor-loop with a regular for-loop.
    pw = PoolWaitbar(MCS, 'Simulation in progress, please wait ...');
    parfor ii = 1:MCS
        [POS_E{ii,1},CPU{ii,1},N_EFF{ii,1}] = main(ii,false,L,J,N);
        increment(pw)
    end
    delete(pw)
    
    % compute performance metrics
    fprintf('N:%d, pos=%.2f [m], Neff=%.2f [%%]\n', ...
        N, ...
        sqrt(mean(cell2mat(POS_E).^2,'all','omitnan')), ...
        mean(cell2mat(N_EFF),'all','omitnan')*100);
    
    % save MCS results
    save(sprintf(filename,J,L,N),'POS_E','CPU','N_EFF')
end

function compute_performance_metrics   
    % create filename
    if strcmp(filesep,'\')
        filename = strjoin([{'results\'} 'PHD_J%dL%dN%d.mat'], filesep);
    else
        filename = strjoin([{'results'} 'PHD_J%dL%dN%d.mat'], filesep);
    end
    
    % define variables
    N = 1;
    L = 5;
    J = 50;
    
    try
        % load data
        load(sprintf(filename,J,L,N),'POS_E','N_EFF');
        
        % compute performance metrics
        fprintf('N:%d, pos=%.2f [m], Neff=%.2f [%%]\n', ...
            N, ...
            sqrt(mean(cell2mat(POS_E).^2,'all','omitnan')), ...
            mean(cell2mat(N_EFF),'all','omitnan')*100);
    catch
        fprintf('Can''t load file, incorrect filename!\n')
    end
end