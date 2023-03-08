function [WPI,MPI,PPI] = gmiplid(obj,col4row,y,P_D,params)
%     This function computes the partitioned GM-ID approximation for J
%     Gaussian mixture components
% 
%     Input:
%        obj     - struct that represent particle n of the PHD-SLAM density
%        col4row - (n_k x J) assignment matrix 
%        y       - {2 x m_k} measurement matrix
%        P_D     - (1 x n_k) detection probabilities
%        params  - simulation parameters
%     
%     Output:
%        WPI      - ((J+1) x 1) GMM weights
%        MPI      - (3 x (J+1)) GMM means
%        PPI      - (3 x 3 x (J+1)) GMM means
%     
%     Author   : Ossi Kaltiokallio
%                Tampere University, Department of Electronics and
%                Communications Engineering
%                Korkeakoulunkatu 1, 33720 Tampere
%                ossi.kaltiokallio@tuni.fi
%     Last Rev : 26/8/2022
%     Tested   : '9.8.0.1359463 (R2020a) Update 1'
%     
%     Copyright notice: You are free to modify, extend and distribute 
%        this code granted that the author of the original code is 
%        mentioned as the original author of the code.  

    % get input variables
    msp = obj.xn;
    Psp = obj.Pn;
    xl = obj.xl;
    Pl = obj.Pl;
    eta = obj.eta;
    
    % determine constant parameters
    R = params.R;
    L = params.L;
    epsilon = params.epsilon;
    lambda_c = params.lambda_c;
    
    % determine dimensions
    J = size(col4row,2);
    n_k = size(xl,2);
    m_k = size(y,2);

    % allocate output variables
    WPI = zeros(1,J+1);
    MPI = zeros(params.xn_dim,J+1);
    PPI = zeros(params.xn_dim,params.xn_dim,J+1);
    
    % compute IPL-OID approximation for each GMM component
    for j = 1:J
        % obtain gamma threshold from a lookup table
        mpi = msp;
        Ppi = Psp;
        
        % Iterations
        l = 0;
        done = ((l >= L) || false);
        while ~done
            % Update iteration counter
            l = l + 1;
            
            % initialize the priors
            ms = msp;
            Ps = Psp;
            
            % compute inverse of linearization covariance
            invPpi = Ppi\eye(size(Ppi));
            
            % compute partitioned OID approximation for l-th iteration
            for i = 1:n_k
                if col4row(i,j) <= m_k
                    % compute moments
                    [my,Hl,Hn] = h_func(mpi,xl(:,i));
                    Py  = Hn*Ppi*Hn' + Hl*Pl(:,:,i)*Hl' + R;
                    Pys = Hn*Ppi;
                    
                    % Linearization
                    A = Pys*invPpi;
                    b = my - A*mpi;
                    Omega = Py - A*Ppi*A';
                    
                    % handle wrap around
                    nu = y(:,col4row(i,j)) - A*ms - b;
                    nu(2) = mod(nu(2) + pi,2*pi) - pi;
                    
                    % update
                    S = A*Ps*A' + Omega;
                    invS = S\eye(size(S));
                    K = Ps*A'*invS;
                    ms = ms + K*nu;
                    Ps = Ps - K*S*K';
                end
            end

            % Check if posterior update was successful, if not, exit loop
            % and use the previous best approximation
            [~, nd] = chol(Ps, 'lower');
            if  nd 
                done = true;
            else
                % compute change in KL divergence
                if l > 1
                    dkls = (trace(Ps\Ppi) - log(det(Ppi)) + log(det(Ps)) - size(ms,1) + (ms - mpi)'/Ps*(ms - mpi))/2;
                    done = dkls < epsilon;
                end
                
                % Update linearization density
                mpi = ms;
                Ppi = Ps;
            end

            % Convergence criteria and tests
            done = (l >= L) || done;
        end
        
        % IPL-OID approximation complete, compute weight for j-th component
        wpi = 0;
        
        for i = 1:n_k
            if col4row(i,j) <= m_k
                % compute moments
                [my,Hl,Hn] = h_func(mpi,xl(:,i));
                Py =  Hn*Ppi*Hn' + Hl*Pl(:,:,i)*Hl' + R;
                
                % handle wrap around
                nu = y(:,col4row(i,j)) - my;
                nu(2) = mod(nu(2) + pi,2*pi) - pi;
                
                % update log-weight with detection
                wpi = wpi + P_D(i) + eta(i) - log(lambda_c) - ...
                    0.5*(params.h_dim*log(2*pi) + log(det(Py)) + nu'*(Py\nu));
            else
                % update log-weight with misdetection
                wpi = wpi + log(1 - exp(P_D(i))) + eta(i);
            end
        end

        
        % append GMM parameters to output
        WPI(1,j) = wpi;
        MPI(:,j) = mpi;
        PPI(:,:,j) = Ppi;
    end
    
    wpi = 0;
    for i = 1:n_k
        wpi = wpi + log(1 - exp(P_D(i))) + eta(i);
    end
    
    WPI(1,end) = wpi;
    MPI(:,end) = msp;
    PPI(:,:,end) = Psp;
end