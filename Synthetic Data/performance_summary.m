function [sim,position,heading,d_gospa] = performance_summary(params,sim)
    % This function calculates the performance metrics

    % Input:
    %    params     - simulation parameters
    %    sim        - struct containing the simulation data
    %
    % Output:
    %    sim        - struct containing the simulation data
    %    position   - (1 x K) position error
    %    heading    - (1 x K) heading error
    %    d_gospa    - (1 x K) GOSPA
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
    
    % if resample flag is false, get trajectory with maximum weight and use
    % that as the estimate
    if ~params.resample
        [~,j] = max(sim.WW(:,end));
        sim.MM = squeeze(sim.XX(:,j,:));
    end

    T = params.T;
    map = sim.map;
    state = sim.state;
    
    % compute position and heading errors
    position = sqrt((state(1,:) - sim.MM(1,:)).^2 + (state(2,:) - sim.MM(2,:)).^2);
    heading = state(3,:) - sim.MM(3,:);
    
    % compute GOSPA 
    d_gospa = nan(1,T);
    if 1
        % only the final map
        [d_gospa(1,T), ~, ~] = GOSPA(map, sim.MM_map{1,T}, 2, 20, 2);
    else
        % all time steps
        FoV = zeros(1,size(map,2));
        for k = 1:T
            
            dx = map(1,:) - state(1,k);
            dy = map(2,:) - state(2,k);
            r_fov = dx.^2 + dy.^2 <= params.fov_range^2;
            b_fov = abs(mod(atan2(dy,dx)-state(3,k) + pi,2*pi) - pi) <= params.fov_angle;
            
            if mod(k,1) == 0
                FoV = FoV | (r_fov & b_fov);
                if k>2 && isempty(sim.MM_map{1,k})
                    d_gospa(1,k) = sum(FoV)*(20^2)/2;
                else
                    [d_gospa(1,k), ~, ~] = GOSPA(map(:,FoV), sim.MM_map{1,k}, 2, 20, 2);
                end
            end
        end
    end
    
    if params.f_mode
        % print performance metrics
        fprintf('N:%d, L:%d, J:%d, pos=%.2f [m], theta=%.2f [deg], gospa=%.2f [m], cpu=%.2f [ms], Time=%.2f [s], Neff=%.2f [%%]\n', ...
            params.N_particle, ...
            params.L, ...
            params.J, ...
            sqrt(mean(position.^2,'omitnan')), ...
            sqrt(mean(heading.^2,'omitnan'))*180/pi, ...
            mean(d_gospa(1,T)), ...
            mean(sim.cpu_time)*1000, ...
            sum(sim.cpu_time), ...
            mean(sim.Neff)*100);

        % illustrate estimated trajectory and map
        figure(1); clf; box on; hold on; grid on;
        plot(map(1,:),map(2,:),'ks','linewidth',1,'markersize',5,'MarkerFaceColor','k');
        plot(sim.MM_map{1,T}(1,:),sim.MM_map{1,T}(2,:),'rs','linewidth',1,'markersize',10);
        plot(state(1,:),state(2,:),'k','linewidth',2)
        plot(sim.MM(1,:),sim.MM(2,:),'r','linewidth',2)
        set(gca,'xlim',[200 1200],'ylim',[200 1200],'XLimMode','manual','YLimMode','manual')
        set(gca,'ticklabelinterpreter','latex','fontsize',16)
        xlabel('x / m','interpreter','latex','fontsize',16)
        ylabel('y / m','interpreter','latex','fontsize',16)
        drawnow

        clear plot_estimate
    end
end