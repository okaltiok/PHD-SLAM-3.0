function [mu,Fx,Fu] = f_func(m,u,dt)
    % This function calculates the state transition and Jacobians of the 
    % vehicle

    % Input:
    %    m      - a (3 x 1) state vector at time k - 1
    %    u      - a (2 x 1) control input at time k - 1
    %    dt     - sampling interval
    %
    % Output:
    %    mu     - (3 x 1) state at time k
    %    Fx     - (3 x 3) Jacobian w.r.t. to vehicle state
    %    Fu     - (3 x 2) Jacobian w.r.t. to control state
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 29/8/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.


    % compute help variables
    s = sin(m(3));
    ss = sin(m(3) + u(2)*dt);
    c = cos(m(3));
    cc = cos(m(3) + u(2)*dt);
    
    % compute mean at time k
    mu = [...
        m(1) - u(1)/u(2)*s + u(1)/u(2)*ss; ...
        m(2) + u(1)/u(2)*c - u(1)/u(2)*cc; ...
        m(3) + u(2)*dt ...
        ];
    
    % compute Jacobian w.r.t. vehicle state
    Fx = [...
        1 0 -u(1)/u(2)*c + u(1)/u(2)*cc; ...
        0 1 -u(1)/u(2)*s + u(1)/u(2)*ss; ...
        0 0 1 ...
        ];
    
    % compute Jacobian w.r.t. control state
    Fu = [-1/u(2)*s + 1/u(2)*ss  u(1)/u(2)^2*s - u(1)/u(2)^2*ss + u(1)/u(2)*dt*cc;
           1/u(2)*c - 1/u(2)*cc -u(1)/u(2)^2*c + u(1)/u(2)^2*cc + u(1)/u(2)*dt*ss;
           0 dt];
end