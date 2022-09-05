function obj = sample(obj,wp,mp,Pp,r_uniform,r_gaussian)
    % This function samples from the Gaussian mixture sampling distribution

    % Input:
    %    obj        - struct that represent particle n of the PHD-SLAM density
    %    wp         - (J x 1) GMM weights
    %    mp         - (3 x J) GMM means
    %    Pp         - (3 x 3 x J) GMM covariances
    %    r_uniform  - uniform random variable
    %    r_gaussian - 3 x 1 normal random variable
    %
    % Output:    
    %    obj        - struct with updated vehicle mean and particle weight  
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
    
    % get number of GMM components
    [xn_dim,J] = size(mp);

    if J > 1
        % normalize log weights
        max_wp = -Inf;
        for i = 1:length(wp)
            if wp(i) > max_wp
                max_wp = wp(i);
            end
        end
        sum_wp = max_wp + log(sum(exp(wp-max_wp)));
        wp = wp-sum_wp;

        % sample categorial distribution to obtain index j
        WP = 0;
        j_idx = 1;
        for j = 1:J
            WP = WP + exp(wp(j));
            if WP >= r_uniform
                j_idx = j;
                break
            end
        end
    else
        wp = 0;
        j_idx = 1;
    end
    
    % check that covariance is positive semidefinite
    [sqrtP, nd] = chol(Pp(:,:,j_idx), 'lower');
    if nd
        % covariance isn't positive semidefinite, do not sample and return
        obj.Pn = zeros(xn_dim);
        obj.iparent = 0;
        return;
    end
    
    % sample from the OID
    xj = mp(:,j_idx) + sqrtP*r_gaussian;

    % compute likelihood of sample w.r.t. to prior 
    d = xn_dim * log(2*pi);
    nu_prior = xj-obj.xn;
    prior = -0.5*(nu_prior'/obj.Pn*nu_prior + d + log(det((obj.Pn))));

    % compute likelihood of sample w.r.t. to proposal
    log_likelihood = zeros(J,1);
    for j = 1:J
        nu_prop = xj-mp(:,j);
        log_likelihood(j) = wp(j) -0.5*(nu_prop'/Pp(:,:,j)*nu_prop + d + log(det((Pp(:,:,j)))));
    end
    prop = log(sum(exp(log_likelihood)));
    
    % update particle weight
    obj.w = obj.w + prior - prop;

    % update UE state and covariance
    obj.xn = xj;
    obj.Pn = zeros(xn_dim);
    obj.iparent = 0;
end