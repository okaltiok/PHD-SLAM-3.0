function [h,Hf,Hv] = h_func(x,lm)
    % This function calculates the mean and Jacobians of the measurements

    % Input:
    %    x      - 3x1 vehicle state
    %    lm     - 2X1 landmark state
    %
    % Output:
    %    h     - mean at time k
    %    Hv    - Jacobian w.r.t. to vehicle 
    %    Hf    - Jacobian w.r.t. to landmark 
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

    % help variables
    dx = lm(1) - x(1);
    dy = lm(2) - x(2);
    d2 = dx^2 + dy^2;
    d  = sqrt(d2);

    % Compute observation
    h = [d; atan2(dy,dx)-x(3)];
    h(2) = mod(h(2) + pi,2*pi) - pi;

    % Jacobian w.r.t. feature states
    if nargout > 1
        Hf = [ dx/d,   dy/d;
            -dy/d2,  dx/d2];
    end

    % Jacobian w.r.t. vehicle state
    if nargout > 2
        Hv = [-dx/d,  -dy/d,   0;
            dy/d2, -dx/d2, -1];
    end
end