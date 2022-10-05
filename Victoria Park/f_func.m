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
    % Last Rev : 26/8/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.

    % Vehicle Configuration
    L = 2.83; %m
    H = 0.75; %m
    a = 3.78; %m
    b = 0.5;  %m

    % compute help variables
    v_c = u(1)/(1-H/L*tan(u(2)));
    sc = a*sin(m(3)) + b*cos(m(3));
    cs = a*cos(m(3)) - b*sin(m(3));
    tmp1 = cos(m(3)) - sc*tan(u(2))/L;
    tmp2 = sin(m(3)) + cs*tan(u(2))/L;
    tmp3 = dt*v_c/(cos(u(2))^2*L);

    % compute mean at time k
    mu = [...
          m(1) + dt*v_c*tmp1; ...
          m(2) + dt*v_c*tmp2; ...
          m(3) + dt*v_c*tan(u(2))/L ...
         ];

    % compute Jacobians
    Fx = [...
          1 0 -dt*v_c*tmp2; ...
          0 1  dt*v_c*tmp1; ...
          0 0 1 ...
         ];

     Fu = [...
         dt*tmp1 -sc*tmp3; ...
         dt*tmp2  cs*tmp3; ...
         dt*tan(u(2))/L tmp3; ...
         ];
end