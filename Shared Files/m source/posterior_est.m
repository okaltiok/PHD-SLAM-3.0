function est = posterior_est(obj,params)
    % This function computes the vehicle and map estimates

    % Input:
    %    obj     - struct that representthe PHD-SLAM density
    %    params  - simulation parameters
    %
    % Output:
    %    est      - struct that contains the estimates
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 26/8/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code. 


    % get index with highest weight
    [~,j] = max([obj.w]);

    % vehicle estimates
    est.m_hat = obj(j).xn;
    est.P_hat = obj(j).Pn;

    % get index of landmarks with weight above a threshold
    idx = obj(j).eta > params.etaT;
    idx_p = obj(j).eta_p > params.etaT;
    nk = sum(idx);
    nk_p = sum(idx_p);
    
    % initialize arrays
    est.mu_hat = zeros(params.xl_dim, nk+nk_p );
    est.C_hat = zeros(params.xl_dim, params.xl_dim, nk+nk_p );
    est.active = cat(2,ones(1,nk),zeros(1,nk_p));

    % append mean and covariance 
    if nk > 0
        est.mu_hat(:,1:nk) =  obj(j).xl(:,idx);
        est.C_hat(:,:,1:nk) =  obj(j).Pl(:,:,idx);
    end

    if nk_p > 0
        est.mu_hat(:,nk+1:end) =  obj(j).xl_p(:,idx_p);
        est.C_hat(:,:,nk+1:end) =  obj(j).Pl_p(:,:,idx_p);
    end
end