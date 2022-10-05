function obj = update(obj,y,params)
    % This function computes the update step of PHD-SLAM  filter

    % Input:
    %    obj        - struct that represent particle n of the PHD-SLAM density at time k | k-1
    %    y          - {2 x m_k} matrix containing the measurements
    %    params     - simulation parameters
    %
    % Output:
    %    obj        - struct that represent particle n of the PHD-SLAM density at time k | k
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
    
    % get number of measurements and return if empty
    m_k = size(y,2);
    if m_k == 0
        obj.birth_y = [];
        return
    end

    % obtain input variables
    xn = obj.xn;
    eta = obj.eta;
    eta_threshold = obj.eta_threshold;
    xl = obj.xl;
    Pl = obj.Pl;

    % number of landmarks and measurements
    n_k = size(xl,2);
    m_k = size(y,2);

    % preallocate output matrices
    hh = zeros(params.h_dim,n_k);
    KK = zeros(params.xl_dim,params.h_dim,n_k);
    PP = zeros(params.xl_dim,params.xl_dim,n_k);
    L = -inf(n_k,m_k);

    % compute detection probabilities
    P_D = adaptiveDetectionProbability(xl, xn, params.fov_range, params.fov_angle, params.P_D);
    
    % compute parameters used in the PHD update
    for i = 1:n_k
        if P_D(i) > 0
            % compute predicted measurement and Jacobian 
            [h,Hl,~] = h_func(xn,xl(:,i));

            % compute innovation covariance and its inverse
            S = Hl*Pl(:,:,i)*Hl' + params.R;
            invS = S\eye(size(S));
            
            % compute log-likelihood coefficient
            c = log(P_D(i)) + eta(i) ...
                - 0.5*(params.h_dim * log(2*pi) + log(det(S)));
                
            % loop through measurements and compute log-likelihoods
            for j = 1:m_k
                nu = y(:,j) - h;
                nu(2) = mod(nu(2) + pi,2*pi) - pi;
                
                % ellipsoidal gating to reduce computational overhead
                d = nu'*invS*nu;
                if d <= params.gating_size
                    L(i,j) = c - 0.5*d;
                end
            end

            % compute Kalman filtering parameters for PHD update
            K =  Pl(:,:,i)*Hl'*invS;
            P = Pl(:,:,i) - K*S*K';
            
            hh(:,i) = h;
            KK(:,:,i) = K;
            PP(:,:,i) = P;
        end
    end


    DA = sum(~isinf(L),1);  % number of landmarks associated per measurement
    NDA = sum(DA,'all');    % number of total DAs

    % reallocate size of output variables
    eta = cat(2,eta,zeros(1,NDA));
    eta_threshold = cat(2,eta_threshold,zeros(1,NDA));
    xl = cat(2,xl,zeros(2,NDA));
    Pl = cat(3,Pl,zeros(2,2,NDA));

    % update weight of misdetection
    for j = 1:n_k
        eta(j) = log(1-P_D(j)) + eta(j);
        eta_threshold(j) = (1-P_D(j)) * eta_threshold(j);
    end

    % compute sum of log weights for each landmark
    eta_tilde_sum = sum(exp(L),1);

    % loop through the measurements and landmarks
    k = n_k+1;
    for i = 1:m_k
        for j = 1:n_k
            % detection
            if ~isinf(L(j,i))
                nu = y(:,i) - hh(:,j);
                nu(2) = mod(nu(2) + pi,2*pi) - pi;

                % update mean and covariance of landmark
                xl(:,k) = xl(:,j) + KK(:,:,j)*nu;
                Pl(:,:,k) = PP(:,:,j);

                % update weight of landmark
                eta(k) = L(j,i) - log(params.lambda_c + eta_tilde_sum(i));
                eta_threshold(k) = params.eta_threshold;

                k = k + 1;
            end
        end
    end

    % compute particle weight
    obj.w = obj.w + sum(log(params.lambda_c + eta_tilde_sum),2);

    % store measurements without a DA that will be used to create
    % birth landmarks
    obj.birth_y = y(:,DA == 0);

    obj.eta = eta;
    obj.eta_threshold = eta_threshold;
    obj.xl = xl;
    obj.Pl = Pl;
end