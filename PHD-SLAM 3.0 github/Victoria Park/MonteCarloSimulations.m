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
    option = 2;
    
    % number of simulations
    N = 100;
    
    switch option
        case 1
            generate_data(N)
        case 2
            perform_mcs(N)
        case 3
            compute_performance_metrics
        case 4
            vp_computational_overhead
        case 5
            vp_parameters_mcs(N)
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

    N_set = 1;%logspace(0,2,3);
    L_set = [0 5 5];
    J_set = [0 1 50];
    resample_flag = [1 0 0];
    
    for n = N_set
        for i = [1 2 3]
            if i == 1
                neff_T = n/2;
            else
                neff_T = 0;
            end
            POS_E = cell(N,1);
            CPU = cell(N,1);
            N_EFF = cell(N,1);
                        
            pw = PoolWaitbar(N, 'Simulation in progress, please wait ...');
            parfor ii = 1:N
                [pos_e,cpu_time,Neff] = main(ii,false,L_set(i),J_set(i),n,resample_flag(i),neff_T);
                               
                POS_E{ii,1} = pos_e;
                CPU{ii,1} = cpu_time;
                N_EFF{ii,1} = Neff;
                
                increment(pw)
            end
            
            delete(pw)
            
            pos_e = cell2mat(POS_E);
            cpu = cell2mat(CPU);
            neff = cell2mat(N_EFF);
            
            fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
                J_set(i), ...
                L_set(i), ...
                n, ...
                sqrt(mean(pos_e.^2,'all','omitnan')), ...
                mean(cpu,'all')*1000, ...
                mean(sum(cpu,2)), ...
                mean(neff,'all')*100);
            
            
            clearvars -except N_eff N POS_E THETA_E GOSPA_D CPU N_EFF J_set L_set N_set resample_flag i ii n
            
            filename = strcat('results/PHD_J',num2str(J_set(i)),'L',num2str(L_set(i)),'N',num2str(n),'.mat');
            save(filename,'POS_E','CPU','N_EFF')
        end
    end
end

function compute_performance_metrics
    % define variables
    N_set = logspace(0,2,3);
    L_set = [0 0 5 5];
    J_set = [0 1 1 50];
    
    % create filename
    filename = 'results/journal results/';
    file = [filename,'PHD_J%dL%dN%d.mat'];

    for n = N_set
        fprintf('\n')
        for i = 1:size(J_set,2)
            % load data
            try
                load(sprintf(file,J_set(i),L_set(i),n));
            catch
                continue
            end
            pos = cell2mat(POS_E);
            cpu = cell2mat(CPU);
            neff = cell2mat(N_EFF);
            
            % compute performance metrics
            fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
                J_set(i), ...
                L_set(i), ...
                n, ...
                sqrt(mean(pos.^2,'all','omitnan')), ...
                mean(cpu,'all')*1000, ...
                mean(sum(cpu,2)), ...
                mean(neff,'all','omitnan')*100);
        end
    end
end


function vp_computational_overhead
    N = 5;
    N_set = fliplr(logspace(0,2,3));
    L_set = [0 0 5 5];
    J_set = [0 1 1 50];
    resample_flag = 0;
    
    for n = N_set
        for i = [4]
            neff_T = 0;
            POS_E = cell(N,1);
            CPU = cell(N,1);
            N_EFF = cell(N,1);
            
            pw = PoolWaitbar(N, 'Simulation in progress, please wait ...');
            for ii = 1:N
                [pos_e,cpu_time,Neff] = main(ii,false,L_set(i),J_set(i),n,resample_flag,neff_T);
                
                POS_E{ii,1} = pos_e;
                CPU{ii,1} = cpu_time;
                N_EFF{ii,1} = Neff;
                
                increment(pw)
%                 pause(max([10 sum(cpu_time,'all')]));
            end
            delete(pw)
            
            pos_e = cell2mat(POS_E);
            cpu = cell2mat(CPU);
            neff = cell2mat(N_EFF);
            
            fprintf('J:%d, L:%d, N:%d, pos=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
                J_set(i), ...
                L_set(i), ...
                n, ...
                sqrt(mean(pos_e.^2,'all','omitnan')), ...
                mean(cpu,'all')*1000, ...
                mean(sum(cpu,2)), ...
                mean(neff,'all','omitnan')*100);
            
            
            clearvars -except N_eff N POS_E THETA_E GOSPA_D CPU N_EFF J_set L_set N_set resample_flag i ii n
            
            filename = strcat('results/PHD_J',num2str(J_set(i)),'L',num2str(L_set(i)),'N',num2str(n),'.mat');
            save(filename,'POS_E','CPU','N_EFF')
        end
    end
end


function vp_parameters_mcs(N)
    n = 10;
    j = 1;
    L_set = 1:1:10;
    for i = L_set
        POS_E = cell(N,1);
        CPU = cell(N,1);
        N_EFF = cell(N,1);
        
        pw = PoolWaitbar(N, 'Simulation in progress, please wait ...');
        parfor ii = 1:N
            [pos_e,cpu_time,Neff] = main(ii,false,i,j,n,false,0);
            POS_E{ii,1} = pos_e;
            CPU{ii,1} = cpu_time;
            N_EFF{ii,1} = Neff;
            increment(pw)
        end
        
        delete(pw)
        
        pos_e = cell2mat(POS_E);
        cpu = cell2mat(CPU);
        neff = cell2mat(N_EFF);
        
        fprintf('J:%d, L:%d, N:%d, pos=%.3f [m], cpu=%.2f [ms], Neff=%.2f [%%]\n', ...
            j, ...
            i, ...
            n, ...
            sqrt(mean(pos_e(:).^2,'omitnan')), ...
            mean(cpu,'all')*1000, ...
            mean(neff,'all')*100);
        
        clearvars -except N POS_E CPU N_EFF L_set i j n
        
        filename = strcat('results/PHD_J',num2str(j),'L',num2str(i),'N',num2str(n),'.mat');
        save(filename,'POS_E','CPU','N_EFF')
    end
end