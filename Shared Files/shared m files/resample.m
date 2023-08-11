function [Neff,flag,obj] = resample(obj, params)
    %RESAMPSTR Stratified resampling
    %
    %   In stratified resampling indices are sampled using random
    %   numbers u_j~U[(j-1)/n,j/n], where n is length of P. Compare
    %   this to simple random resampling where u_j~U[0,1]. See,
    %   Kitagawa, G., Monte Carlo Filter and Smoother for
    %   Non-Gaussian Nonlinear State Space Models, Journal of
    %   Computational and Graphical Statistics, 5(1):1-25, 1996. 
    %
    % Copyright (c) 2003-2004 Aki Vehtari
    %
    % This software is distributed under the GNU General Public 
    % Licence (version 2 or later); please refer to the file 
    % Licence.txt, included with the software, for details.
    %
    % The algorithm below has been modified for the purposes of PHD-SLAM
    %
    % Input:
    %    obj    - a (1 x N) struct that represent the PHD-SLAM density
    %    params - simulation parameters
    %    
    % Output:
    %    Neff   - Effective sample size
    %    flag   - true if resampling is performed
    %    obj    - a resampled (1 x N) struct that represent the PHD-SLAM density
    %
    % Modified by : Ossi Kaltiokallio
    %               Tampere University, Department of Electronics and
    %               Communications Engineering
    %               Korkeakoulunkatu 1, 33720 Tampere
    %               ossi.kaltiokallio@tuni.fi
    % Last Rev    : 1/3/2022
    % Tested      : Matlab version 9.8.0.1359463 (R2020a) Update 1
    
    
    N = params.N_particle;
    flag = 0;
    if N == 1
        Neff = 1;
        return
    end
    
    % normalize log weights 
    log_w = [obj.w];
    max_log_w = -Inf;
    for i = 1:N
        if log_w(i) > max_log_w
            max_log_w = log_w(i);
        end
    end
    log_sum_w = max_log_w + log(sum(exp(log_w-max_log_w)));
    log_w = log_w-log_sum_w;
    
    % convert to linear scale
    w = exp(log_w);

	% compute effective number of samples
    Neff = 1/sum(w.^2);
    if params.resample && Neff <= params.N_eff
        flag = 1;
        wn = w.*N;
        
        s=zeros(N,1);
        r=rand(N,1);
        k=0;
        c=0;
        for i=1:N
            c=c+wn(i);
            if c>=1
                a=floor(c);
                c=c-a;
                s(k+(1:a))=i;
                k=k+a;
            end
            if k<N && c>=r(k+1)
                c=c-1;
                k=k+1;
                s(k)=i;
            end
        end

        % replace particles with low weight
        obj_temp = obj;
        log_w = log(1/N);
        for i = 1:N
            obj(i) = obj_temp(s(i));
            obj(i).w = log_w;
            obj(i).iparent = s(i);
        end
    end
end