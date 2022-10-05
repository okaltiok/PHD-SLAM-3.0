function [obj,est] = rbpf_phd(obj,y,u,t,params)
    % This function performs one iteration of the PHD-SLAM filter recursion 
    % including prediction, computing the proposal, sampling, update, map
    % pruning, estimation and resampling

    % Input:
    %    obj    - a (1 x N) struct that represent the PHD-SLAM density at time k - 1
    %    y      - a (2 x m_k) matrix that contains the received measurements
    %    u      - a (2 x t_k) control input 
    %    t      - a (1 x t_k) time vector 
    %    params - simulation parameters
    %
    % Output:
    %    obj   - updated PHD-SLAM density at time k
    %    est   - a struct that contains the vehicle and map estimates
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
    
    tic

    N = params.N_particle;
    iPrevParent = 0;
    
    % generate random numbers used in sampling
    r_uniform = rand(1,N);
    r_gaussian = randn(params.xn_dim,N);

    % for each particle do
    for i = 1:N
        % predict and make a temporary copy of the particle
        obj_i = predict(obj(i),u,t,params);
        
        % compute parameters of the proposal
        iParent = obj_i.iparent;
        if iParent == 0 || iParent ~= iPrevParent
            if params.J == 0 || size(obj_i.xl,2) == 0 || size(y,2) == 0
                % when there are no features nor measurements, use the
                % motion model as the sampling distribution
                wp = 0;
                mp = obj_i.xn;
                Pp = obj_i.Pn;
            else              
                % compute cost matrix
                [P_D, L] = cost_matrix(obj_i, y, params);
                
                % solve ranked assignment problem
                col4rowBest = kBest2DAssign(L,params.J,true,params.DA_threshold);

                % compute partitioned GM-OID approximation
                [wp,mp,Pp] = mhrbiploid(obj_i,double(col4rowBest),y,log(P_D),params);
            end
            iPrevParent = iParent;
        end
        
        % sample from proposal
        obj_i = sample(obj_i,wp,mp,Pp,r_uniform(i),r_gaussian(:,i));

        % update PHD parameters
        obj_i = update(obj_i,y,params);

        % prune and merge PHD components, and store the particle in the 
        % PHD-SLAM density object
        obj(i) = hypothesisReductionAlgorithm(obj_i,params);
    end

    % estimate and resample
    est = posterior_est(obj,params);
    [Neff,obj] = resample(obj, params);

    est.Neff = Neff/N;
    est.dt = toc;
end