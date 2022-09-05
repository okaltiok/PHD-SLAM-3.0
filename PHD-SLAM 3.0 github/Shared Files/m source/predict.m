function obj = predict(obj,u,t,params)
    % This function implements the PHD prediction step including the birth 
    % process that creates new landmarks the propagation of the vehicle
    % according to the kinematic model and control inputs

    % Input:
    %    obj        - struct that represent particle n of the PHD-SLAM density at time k-1
    %    u          - {2 x T} control matrix
    %    t          - {1 x T} time vector
    %    params     - simulation parameters
    %
    % Output:
    %    obj        - struct that represent particle n of the PHD-SLAM density at time k | k-1
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
    
    % get input variables
    xn = obj.xn;
    Pn = obj.Pn;
    xl = obj.xl;
    Pl = obj.Pl;
    eta = obj.eta;
    eta_threshold = obj.eta_threshold;


    % create new landmarks
    m_k = size(obj.birth_y,2);
    if m_k > 0
        % allocate matrices for new landmarks
        xb = zeros(params.xl_dim,m_k);
        Pb = zeros(params.xl_dim,params.xl_dim,m_k);
        R = params.R;
        
        % create new landmark for every measurement
        for j = 1:m_k
            % help variables
            r = obj.birth_y(1,j);
            b = obj.birth_y(2,j);
            s= sin(xn(3) + b);
            c= cos(xn(3) + b);
            
            % Jacobian
            Gz= [c -r*s;
                 s  r*c];
            
            % mean
            xb(:,j) = [xn(1) + r*c;
                xn(2) + r*s];
           
            % covariance
            Pb(:,:,j) = Gz*R*Gz';
        end

        % concatenate old landmarks with new
        xl = cat(2,xl,xb);
        Pl = cat(3,Pl,Pb);
        eta = cat(2,eta,log(params.P_B)*ones(1,m_k));
        eta_threshold = cat(2,eta_threshold,params.eta_threshold*ones(1,m_k));
    end

    %  propagate vehicle
    T = size(t,2);
    if T > 1
        Qu = params.Qu;
        Pn = Pn + params.Qn;
        for k = 2:T
            dt = t(k) - t(k-1);
            [xn,Fx,Fu] = f_func(xn,u(:,k-1),dt);
            Pn = Fx*Pn*Fx' + Fu*Qu*Fu';
        end
    end
    
    obj.xn = xn;
    obj.Pn = Pn;
    obj.xl = xl;
    obj.Pl = Pl;
    obj.eta = eta;
    obj.eta_threshold = eta_threshold;
    obj.birth_y = [];
end