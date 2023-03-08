function [P_D, L] = cost_matrix(obj, Y, params)
    % This function computes the cost matrix and detection probability

    % Input:
    %    obj - struct that represent particle n of the PHD-SLAM density
    %    Y          - {2 x m_k} matrix containing the measurements
    %    params     - simulation parameters
    %
    % Output:
    %    P_D        - (1 x n_k) containing the detection probabilities
    %    L          - (n_k x [n_k + m_k]) cost matrix
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
    
    xn = obj.xn;
    Pn = obj.Pn;
    xl = obj.xl;
    Pl = obj.Pl;
    eta = obj.eta;
    
    % determine dimensions
    n_k = size(xl,2);
    m_k = size(Y,2);
  
    L = -inf(n_k,m_k+n_k);

    % compute detection probabilities
    P_D = adaptiveDetectionProbability(xl, xn, params.fov_range, params.fov_angle, params.P_D);
    
    % loop through landmarks
    for i = 1:n_k       
        L(i,m_k+i) = log(1 - P_D(i)) + eta(i);
        
        % compute predicted measurement and Jacobian 
        if P_D(i) > 0
            [h,Hl,Hn] = h_func(xn,xl(:,i));

            % compute innovation covariance
            S = Hl*Pl(:,:,i)*Hl' + Hn*Pn*Hn' + params.R;
            
            % compute normalized residual
            invS = S\eye(size(S));

            % compute log-likelihood coefficient
            c = log(P_D(i)) - log(params.lambda_c) + eta(i) ...
                - 0.5*(params.h_dim * log(2*pi) + log(det(S)));
            
            % loop through measurements and compute log-likelihoods
            for j = 1:m_k
                % residual
                nu = Y(:,j) - h;
                nu(2) = mod(nu(2) + pi,2*pi) - pi;
                
                d = nu'*invS*nu;
                
                % ellipsoidal gating to reduce computational overhead
                if d <= params.gating_size
                    L(i,j) = c - 0.5*d;
                end
            end
        end
    end
end
