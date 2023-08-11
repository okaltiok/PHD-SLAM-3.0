function obj =  hypothesisReductionAlgorithm(obj,params)
    % This function performs GM component pruning and merging, and
    % determines the components that are active and passive. An active
    % component is within the FOV or close to it, and it is passive 
    % otherwise.

    % Input:
    %    obj        - a struct that represent the updated particle (i) of the PHD-SLAM density at time k
    %    params     - simulation parameters

    % Output:
    %    obj        - obj with pruned and merged GM components
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 30/3/2022
    % Tested   : Matlab version 9.8.0.1359463 (R2020a) Update 1
    %
    % Copyright notice: You are free to modify, extend and distribute
    %    this code granted that the author of the original code is
    %    mentioned as the original author of the code.

    % get parameters of the Gaussian mixture
    eta = obj.eta;
    xl = obj.xl;
    Pl = obj.Pl;
    xl_dim = size(xl,1);
    
    % GM reduction thresholds
    w_min = params.w_min;                               % pruning threshold
    merging_threshold = params.merging_threshold;       % merging threshold
    
    % number of GM components
    n_k = size(eta,2);
    
    % Sort weights in descending order using Bubble sort
    sorted_idx = 1:n_k;
    endofarray = n_k;
    for i = 1:n_k
        % Pass up to the last un-sorted element.
        for j = 1:n_k-i
            % If elements are in the wrong order, swap them.
            if(eta(j+1)>eta(j))
                temp = eta(j);
                eta(j) = eta(j+1);
                eta(j+1) = temp;
                
                temp = sorted_idx(j);
                sorted_idx(j) = sorted_idx(j+1);
                sorted_idx(j+1) = temp;
            end
        end

        if i < n_k && (eta(j+1) < w_min || isinf(eta(j+1)))
            endofarray = j;
        end
    end
    
    % sort the other parameters as well
    xl = xl(:,sorted_idx);
    Pl = Pl(:,:,sorted_idx);
    

    % merge GM components that are close to one another
    merged_idx = zeros(1,endofarray);
    cluster_idx = zeros(1,endofarray);
    el = 1;
    for i = 1:endofarray
        if merged_idx(i)
            continue
        end
        
        % find components that are close to one another
       
        cluster_idx(i) = 1;
        merged_idx(i) = 1;
        n = 1;
        
        invP = Pl(:,:,i)\eye(xl_dim);
        
        for j = i+1:endofarray
            
            nu = xl(:,i) - xl(:,j);
            
            if nu'*invP*nu < merging_threshold
                cluster_idx(j) = 1;
                merged_idx(j) = 1;
                n = n + 1;
            end
        end
        
        % normalize log weights of the components that are merged
        log_w = zeros(1,n);
        idx = zeros(1,n);
        log_sum_w = 0;
        l = 1;
        for j = i:endofarray
            if cluster_idx(j)
                log_w(l) = eta(j);
                idx(l) = j;
                log_sum_w = log_sum_w + exp(log_w(l)-log_w(1));
                l = l + 1;
                cluster_idx(j) = 0;
            end
        end
        log_sum_w = log_w(1) + log(log_sum_w);
        log_w = log_w-log_sum_w;
        

        % merge components using Gaussian moment matching
        if n == 1
            m = xl(:,idx(1));
            P = Pl(:,:,idx(1));
        else
            w = exp(log_w);
            m = zeros(xl_dim,1);
            P = zeros(xl_dim,xl_dim);
            for j = 1:n
                m = m + w(j).*xl(:,idx(j));
            end
            for j = 1:n
                l = idx(j);
                nu = (m - xl(:,l));
                for jj = 1:xl_dim
                    for jjj = 1:xl_dim
                        P(jj,jjj) = P(jj,jjj) + w(j).*Pl(jj,jjj,l) + w(j).*nu(jj)*nu(jjj);
                    end
                end
            end
        end

        
        if trace(P) > params.merging_threshold
            log_sum_w = log(1e-18);
        end
        
        % append the merged components to the parameters
        eta(el) = log_sum_w;
        xl(:,el) = m;
        Pl(:,:,el) = P;
        
        el = el+1;
    end
    
    % determine active and passive components
    xn = obj.xn;
    xlp = obj.xl_p;
    xl_k = el-1;
    xlp_k = size(xlp,2);
    
    FOV_RANGE2 = (params.fov_range + params.range_buffer)^2;
    FOV_ANGLE = params.fov_angle + params.angle_buffer;
    
    % determine active components that are inside/outside FOV
    xl_fov = ones(xl_k,1);
    for i = 1:xl_k
        dx = xl(1,i) - xn(1);
        dy = xl(2,i) - xn(2);
        
        if dx^2 + dy^2 > FOV_RANGE2
            xl_fov(i,1) = 0;
        elseif abs(mod(atan2(dy,dx)-xn(3) + pi,2*pi) - pi) > FOV_ANGLE
            xl_fov(i,1) = 0;
        end
        
        if xl_fov(i,1) == 0 && eta(i) <= w_min
            xl_fov(i,1) = -1;
        end
    end
    
    % determine passive components that are inside/outside FOV
    xlp_fov = ones(xlp_k,1);
    for i = 1:xlp_k
        dx = xlp(1,i) - xn(1);
        dy = xlp(2,i) - xn(2);
        
        if dx^2 + dy^2 > FOV_RANGE2
            xlp_fov(i,1) = 0;
        elseif abs(mod(atan2(dy,dx)-xn(3) + pi,2*pi) - pi) > FOV_ANGLE
            xlp_fov(i,1) = 0;
        end
    end
    
    move_outside_fov = sum(xl_fov==0,1);
    move_inside_fov = sum(xlp_fov==1,1);

    % store active components
    if move_inside_fov == 0
        obj.xl = xl(:,xl_fov==1);
        obj.Pl = Pl(:,:,xl_fov==1);
        obj.eta = eta(:,xl_fov==1);
    else
        obj.xl = cat(2,xl(:,xl_fov==1),xlp(:,xlp_fov==1));
        obj.Pl = cat(3,Pl(:,:,xl_fov==1),obj.Pl_p(:,:,xlp_fov==1));
        obj.eta = cat(2,eta(:,xl_fov==1),obj.eta_p(:,xlp_fov==1));
    end
    
    % store passive components
    if move_outside_fov == 0
        obj.xl_p = xlp(:,xlp_fov==0);
        obj.Pl_p = obj.Pl_p(:,:,xlp_fov==0);
        obj.eta_p = obj.eta_p(:,xlp_fov==0);
    else
        obj.xl_p = cat(2,xlp(:,xlp_fov==0),xl(:,xl_fov==0));
        obj.Pl_p = cat(3,obj.Pl_p(:,:,xlp_fov==0),Pl(:,:,xl_fov==0));
        obj.eta_p = cat(2,obj.eta_p(:,xlp_fov==0),eta(:,xl_fov==0));
    end
    
    if size(obj.xl_p,2) == 0
        obj.xl_p = zeros(2,0);
        obj.Pl_p = zeros(2,2,0);
        obj.eta_p = zeros(1,0);
    end
end