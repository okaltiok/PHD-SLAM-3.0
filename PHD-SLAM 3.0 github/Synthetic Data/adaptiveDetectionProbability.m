function P_D = adaptiveDetectionProbability(xl, xn, FOV_RANGE, FOV_ANGLE, DETECTION_PROBABILITY)
    % This function calculates the detection probability of xl observed from xn

    % Input:
    %    xl                     - a (2 x n_k) matrix that contains the landmark states
    %    xn                     - a (3 x 1) state vector
    %    FOV_RANGE              - maximum scanning range of sensor
    %    FOV_ANGLE              - maximum scanning angle of sensor
    %    DETECTION_PROBABILITY  - detection probability inside FOV
    %
    % Output:
    %    P_D                    - a (1 x n_k) detection probability vector
    %
    % Author   : Ossi Kaltiokallio
    %            Tampere University, Department of Electronics and
    %            Communications Engineering
    %            Korkeakoulunkatu 1, 33720 Tampere
    %            ossi.kaltiokallio@tuni.fi
    % Last Rev : 1/9/2022
    % Tested   : '9.8.0.1359463 (R2020a) Update 1'
    %
    % Copyright notice: You are free to modify, extend and distribute 
    %    this code granted that the author of the original code is 
    %    mentioned as the original author of the code.   
    
    n_k = size(xl,2);
    P_D = zeros(1,n_k);
    
    if n_k > 0
        FOV_RANGE2 = FOV_RANGE^2;
        
        % compute squared range and angles
        dx = xl(1,:) - xn(1);
        dy = xl(2,:) - xn(2);
        range2 = dx.^2 + dy.^2;
        theta = mod(atan2(dy,dx)-xn(3) + pi,2*pi) - pi;

        for i = 1:n_k
            if range2(i) > FOV_RANGE2
                continue
            elseif abs(theta(i)) > FOV_ANGLE
                continue
            else
                % within FOV, set probability of detection
                P_D(i) = DETECTION_PROBABILITY;
            end
        end
    end
end