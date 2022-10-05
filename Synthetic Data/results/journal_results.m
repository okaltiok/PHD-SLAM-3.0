function journal_results
    str = {'PHD 1.0','PHD 2.0b','PHD 3.0'};
    J_set = [0 1 50];
    L_set = [0 5 5];
    N_set = logspace(0,2,3);
    
    filename = 'without resampling/';

    file = [filename,'PHD_J%dL%dN%d.mat'];

    pos_array = zeros(size(J_set,2),size(N_set,2));
    gospa_array = zeros(size(J_set,2),size(N_set,2));
    for n = 1:size(N_set,2)
        for i = 1:size(J_set,2)
            try
                load(sprintf(file,J_set(i),L_set(i),N_set(n)));
            catch
                continue
            end
            pos_e = cell2mat(POS_E);
            theta_e = cell2mat(THETA_E);
            gospa_d = cell2mat(GOSPA_D);
            cpu = cell2mat(CPU);
            neff = cell2mat(N_EFF);
            resamp = (sum(neff <= 0.5,'all')./numel(neff));
            
            pos_array(i,n) = sqrt(mean(pos_e.^2,'all'));
            gospa_array(i,n) = mean(gospa_d(:,end),'all');
            
            fprintf('%s - pos=%.2f (%.2f) [m], theta=%.2f (%.2f) [deg], gospa=%.2f (%.2f)  [m], Neff=%.2f [%%], Resampling=%.2f [%%], cpu=%.2f (%.2f) [ms]. Time=%.2f [a]\n', ...
                str{i}, ...
                sqrt(mean(pos_e.^2,'all')), std(pos_e,0,'all'), ...
                sqrt(mean(theta_e.^2,'all'))*180/pi, std(theta_e,0,'all')*180/pi, ...
                mean(gospa_d(:,end),'all'), std(gospa_d(:,end),0,'all'), ...
                mean(neff,'all')*100, ...
                resamp*100, ...
                mean(cpu,'all')*1000, std(cpu,0,'all')*1000,mean(sum(cpu,2)));
            
        end
    end
    
    fprintf('\n')
    
    for i = 1:size(J_set,2)
        fprintf('%s ',str{i});
        for j = 1:size(N_set,2)
                fprintf('& $%.1f$ & $%.1f$ ',pos_array(i,j),gospa_array(i,j));
            end
        fprintf('\\\\\n')
    end
    
    fprintf('\n')
    
    filename = 'with resampling/';
    file = [filename,'PHD_J%dL%dN%d.mat'];
    for n = 1:size(N_set,2)
        for i = 1:size(J_set,2)
            try
                load(sprintf(file,J_set(i),L_set(i),N_set(n)));
            catch
                continue
            end
            pos_e = cell2mat(POS_E);
            theta_e = cell2mat(THETA_E);
            gospa_d = cell2mat(GOSPA_D);
            
            neff = cell2mat(N_EFF);
            resamp = (sum(neff <= 0.5,'all')./numel(neff));
            
            try
                load(sprintf(['with resampling/computational complexity 2/','PHD_J%dL%dN%d.mat'],J_set(i),L_set(i),N_set(n)));
                cpu = cell2mat(CPU);
            catch
                cpu = cell2mat(CPU);
%                 continue
            end
            
            if 1
                fprintf('%s - pos=%.2f (%.2f) [m], theta=%.2f (%.2f) [deg], gospa=%.2f (%.2f)  [m], Neff=%.2f [%%], Resampling=%.2f [%%], cpu=%.2f (%.2f) [ms]. Time=%.2f [a]\n', ...
                    str{i}, ...
                    sqrt(mean(pos_e.^2,'all')), std(pos_e,0,'all'), ...
                    sqrt(mean(theta_e.^2,'all'))*180/pi, std(theta_e,0,'all')*180/pi, ...
                    mean(gospa_d(:,end),'all'), std(gospa_d(:,end),0,'all'), ...
                    mean(neff,'all')*100, ...
                    resamp*100, ...
                    mean(cpu,'all')*1000, std(cpu,0,'all')*1000,mean(sum(cpu,2)));
            else
                fprintf('%s & $%d$ & $%.2f \\pm %.2f$  & $%.2f \\pm %.2f$ & $%.2f \\pm %.2f$ & $%.2f$ & $%.2f$ & $%.2f \\pm %.2f$ \\\\ \n',...
                    str{i}, ...
                    N_set(n), ...
                    sqrt(mean(pos_e.^2,'all')), std(pos_e,0,'all'), ...
                    sqrt(mean(theta_e.^2,'all'))*180/pi, std(theta_e,0,'all')*180/pi, ...
                    mean(gospa_d(:,end),'all'), std(gospa_d(:,end),0,'all'), ...
                    mean(neff,'all')*100, ...
                    resamp*100, ...
                    mean(cpu,'all')*1000, std(cpu,0,'all')*1000);
                
            end
        end
    end
    
    