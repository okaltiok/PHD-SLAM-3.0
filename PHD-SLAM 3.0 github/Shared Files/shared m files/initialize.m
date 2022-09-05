function obj = initialize(sim,params)
    % Initialize PHD-SLAM density

    % Input:
    %    sim    - struct containing the initial estimates
    %    params - simulation parameters
    %
    % Output:
    %    obj    - a (1 x N) struct that represent the PHD-SLAM density
    %
    %    The fields of particle j are:
    %    xn             - state of the vehicle
    %    Pn             - covariance of the vehicle
    %    w              - particle weight
    %    iParent        - index of the parent particle
    %    eta            - (1 x M) vector of landmark weights
    %    eta_threshold  - (1 x M) vector of landmark weight threshold
    %    xl             - (2 x M) matrix of landmark states
    %    Pl             - (2 x 2 x M) matrix of landmark covariances
    %    eta_p          - (1 x M) vector of landmark weights (outside FOV)
    %    eta_threshold_p  - (1 x M) vector of landmark weight threshold (outside FOV)
    %    xl_p           - (2 x M) matrix of landmark states (outside FOV)
    %    Pl_p           - (2 x 2 x M) matrix of landmark covariances (outside FOV)
    %    birth_y        - (2 x C) matrix containing measurements used to create new landmarks
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
    

    % obj is a (1 x N) struct that represent the PHD-SLAM density
    obj = repmat( ...
        struct('xn',sim.x0, ...
        'Pn',sim.P0, ...
        'w',log(1/params.N_particle), ...
        'iparent',1, ...
        'eta',[], ...
        'eta_threshold',[], ...
        'xl',[], ...
        'Pl',[], ...
        'eta_p',[], ...
        'eta_threshold_p',[], ...
        'xl_p',[], ...
        'Pl_p',[], ...
        'birth_y',[]) ...
        ,1,params.N_particle);
end
