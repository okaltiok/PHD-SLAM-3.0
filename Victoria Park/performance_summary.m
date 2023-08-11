function position = performance_summary(params,sim)
    % This function calculates the performance metrics

    % Input:
    %    params     - simulation parameters
    %    sim        - struct containing the simulation data
    %
    % Output:
    %    sim        - struct containing the simulation data
    %    position   - (1 x K) position error
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
    
    % initialize (1 x K) RMS position error vector
    position = nan(1,params.T);

    % for each time instant compute the RMSE
    for j = 1:params.T
        if ~isempty(sim.o_time{j})
            [val,l] = min(abs(sim.gps_time-sim.o_time{j}));
            
            % only use GPS values for which we have an estimate within
            % 25 milliseconds 
            if val <= 0.025
                position(j) = sqrt(sum((sim.state(:,l) - sim.MM(1:2,j)).^2,1));
            end
        end
    end

    if params.f_mode
        % print performance metrics
        pos_rmse = sqrt(mean(position.^2,'omitnan'));
        fprintf('N:%d, L:%d, J:%d, pos=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%], Resample==%.2f [%%]\n', ...
            params.N_particle, ...
            params.L, ...
            params.J, ...
            pos_rmse, ...
            mean(sim.cpu_time)*1000, ...
            sum(sim.cpu_time), ...
            mean(sim.Neff)*100, ...
            sum(sim.resample)/params.T*100);

        % illustrate estimated trajectory and map
        figure(1); clf; cla; hold on; box on; grid on;
        title('')
        x = sim.MM(1,:);
        y = sim.MM(2,:);
        x(x == 0) = nan;
        y(y == 0) = nan;

        for i = size(sim.MM_map,2):-1:1
            if ~isempty(sim.MM_map{1,i})
                map = sim.MM_map{1,i};
                break
            end
        end

        plot(map(1,:),map(2,:),'b+','linewidth',1,'markersize',8,'visible','on');
        plot(sim.state(1,:),sim.state(2,:),'k.','markersize',8)
        plot(x,y,'r','linewidth',3)
        axis([-210 90 -130 220])
        set(gca,'ticklabelinterpreter','latex','fontsize',16)
        xlabel('x / m','interpreter','latex','fontsize',16)
        ylabel('y / m','interpreter','latex','fontsize',16)
        drawnow

        clear plot_estimate
    end
end