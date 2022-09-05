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

    % get current folder and define path name
    folderParts = strsplit(pwd, filesep);
    meas_path = strjoin([folderParts 'measurements\VP data set'], filesep);
    
    % add measurement path to measurements
    if ~contains(path,meas_path), addpath(meas_path); end

    % create data and save
    for ii = 1:N
        clear('params','model','sim','obj')
        [params,sim,filename]  = simulation_setup(ii,'save');
        save(filename, 'params', 'sim')
        fprintf('round %d\n',ii)
    end
    
    % remove created paths
    rmpath(meas_path);
end

function perform_mcs(N)
    folderParts = strsplit(pwd, filesep);
    shared_m_path = strjoin([folderParts(1:end-1) 'Shared Files\shared m files'], filesep);
    if ~contains(path,shared_m_path), addpath(shared_m_path); end

    L = 5;
    J = 50;
    N_particle = logspace(0,2,3);
    
    for n = N_particle

        POS_E = cell(N,1);
        CPU = cell(N,1);
        N_EFF = cell(N,1);
        
        pw = PoolWaitbar(N, 'Simulation in progress, please wait ...');
        parfor i = 1:N
            [POS_E{i,1},CPU{i,1},N_EFF{i,1}] = main(i,false,L,J,n);
            increment(pw)
        end
            
        delete(pw)

        fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], cpu=%.2f [ms]\n', ...
            J, ...
            L, ...
            n, ...
            sqrt(mean(cell2mat(POS_E).^2,'all','omitnan')), ...
            mean(cell2mat(CPU),'all')*1000);
        
        clearvars -except N POS_E CPU N_EFF J L N_particle i n
        
        filename = strcat('results/PHD_J',num2str(J),'L',num2str(L),'N',num2str(n),'.mat');
        save(filename,'POS_E','CPU','N_EFF')
    
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
        pos = cell2mat(POS_E);
        cpu = cell2mat(CPU);
        
        % compute performance metrics
        fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], cpu=%.2f [ms]\n', ...
            J, ...
            L, ...
            n, ...
            sqrt(mean(pos.^2,'all','omitnan')), ...
            mean(cpu,'all')*1000);

    end
end