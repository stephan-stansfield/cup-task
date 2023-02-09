function [sys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,forwardF,ver,impedance,printSys)

    % DEBUG
%     disp("Inside sysCreate")

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % First create system then get damping ratio & natural frequency info

    % Simulate system including hand impedance
    if impedance

        % Multi-mode model
        if ver == "full"
        
            A = [      0,                0,               1,      0 ;       % * [    x
                       0,                0,               0,      1 ;       %      theta
                    -k / M,          g * m / M,        -b / M,    0 ;       %      x_dot
                 k / (l * M), - g / l * (m + M) / M, b / (l * M), 0 ];      %    theta_dot ]
    
            % Input U(t) includes feedforward force and half of impedance force
            B = [      0 ;                                                  % * U(t)
                       0 ;
                     k / M ;
                  -k / (l * M) ];
    
            C = eye(4);
    
            D = zeros(4, 1);
            
            Arigid = [       0,       0,      1,       0 ;                  % * [    x 
                             0,       0,      0,       1 ;                  %      theta
                        -k / (m + M), 0, -b / (m + M), 0 ;                  %      x_dot    
                             0,       0,      0,       0 ];                 %    theta_dot ]
                        
            Brigid = [    0 ;                                               % * U(t)
                          0 ;
                      k / (m + M);
                          0 ];
        
        % Single mode models
        else 

            Arigid = zeros(4);
            Brigid = zeros(4, 1);

            % Simplification: Lumped mass cup and ball
            if ver == "rigid body"
    
                A = [       0,       0,      1,       0 ;                   % * [    x 
                            0,       0,      0,       1 ;                   %      theta
                       -k / (m + M), 0, -b / (m + M), 0 ;                   %      x_dot    
                            0,       0,      0,       0 ];                  %    theta_dot ]
                        
                B = [    0 ;                                                % * U(t)
                         0 ;
                      k / (m + M) ;
                         0 ];

                C = eye(4);

                D = zeros(4, 1);

            % Simplification: No hand impedance
            elseif ver == "no impedance"

                A=[0,    0,   1, 0 ;                                        % * [    x                
                   0,    0,   0, 1 ;                                        %      theta
                   0,    0,   0, 0 ;                                        %      x_dot    
                   0, -g / l, 0, 0 ];                                       %    theta_dot ]

                B=[   0;                                                    % * ac (shaped x acceleration)
                      0;
                      1;
                    -1/l];

                C = eye(4);
        
                D = zeros(4, 1);
            
            elseif ver == "slow"
                A = [];
                B = [];
                C = eye(4);
                D = zeros(4, 1);
    
            elseif ver == "fast"
                A = [];
                B = [];
                C = eye(4);
                D = zeros(4, 1);
    
            end

        end
    
    % No hand impedance (used in Guang et al, 2019)
    else
        % Linearized 1D pendulum below cart, no stiffness or damping
        A=[0   0  1 0;                                                      % *[x theta x_dot theta_dot]'
           0   0  0 1;
           0   0  0 0;
           0 -g/l 0 0];

        % Single input "u" is desired acceleration
        B=[   0;                                                            % *x''_des'
              0;
              1;
            -1/l];

        % 4 outputs: cart position, ball angle, cart vel, & ball angular
        % velocity. Ball angle converted from radians to degrees.
        C = eye(4);

        D = zeros(4, 1);

        states={'x' 'theta' 'x_{dot}' 'theta_{dot}'};
        
        Arigid =   [0 0 1 0;
                    0 0 0 1;
                    0 0 0 0;
                    0 0 0 0];
                
        Brigid =   [0;
                    0;
                    1;
                    0];

        % Assign natural frequency and damping ratio to output variables
        [w,z] = damp(sys);
        omega = w(3);
        zeta = z(3);
        Td = 2*pi/omega;
        overdamped = false;
        
    end

    states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
    inputs = {'u'};
    outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
    sys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, ...
        'outputname', outputs);
    sysRigid = ss(Arigid, Brigid, C, D, 'statename', states, ...
            'inputname', inputs, 'outputname', outputs);
    
    % Print out system transfer function, poles, and zeros
    if printSys
        H = ss2tf(A, B, C, D)

        disp(['"zeros" contains the system zeros in its columns. Number of ' ...
        'columns in zeros corresponds to number of outputs. "poles" is a ' ...
        'column vector containing the coefficients of the poles of the ' ...
        'poles of the system.']);
        [zs, poles, ~] = ss2zp(A, B, C, D)
        
        % Show pole-zero map
        figure();
        pzmap(sys)

    end

    %%%%%%%
    % Get natural frequency and damping ratio info
    [w,z] = damp(sys);

    % Multi-mode system
    if ver == "full"
        % Assign natural frequency and damping ratio of two modes to variables
        w1 = w(1);
        w2 = w(3);
        z1 = z(1);
        z2 = z(3);

        % Input shaper using LOWER frequency
        zeta1 = z1;
        omega1 = w1;

        if zeta1 >= 1
            % If overdamped, raise flag and set other outputs to zero
            overdamped = true;
            Td1 = 0;
            Td2 = 0;
            zeta1 = 0;
            zeta2 = 0;
        else
            omega_d1 = omega1 * sqrt(1 - zeta1 ^ 2);                    % Damped natural frequency, lower mode (rad/s)
            Td1 = 2 * pi / omega_d1;                                    % Damped period of vibration, lower mode (s)
            overdamped = false;
        end

        % Input shaper using HIGHER frequency
        zeta2 = z2;
        omega2 = w2;

        % If first mode isn't overdamped, check second mode
        if not(overdamped)
            if zeta2 >= 1
                % If overdamped, raise flag and set other outputs to zero
                overdamped = true;
                Td1 = 0;
                Td2 = 0;
                zeta1 = 0;
                zeta2 = 0;
            else
                omega_d2 = omega2 * sqrt(1 - zeta2 ^ 2);                % Damped natural frequency, higher mode (rad/s)
                Td2 = 2 * pi / omega_d2;                                % Damped period of vibration, higher mode (s)
                overdamped = false;
            end
        end
        
        % Set single-mode period & damping ratio outputs to zero
        Td = 0;
        zeta = 0;
        
        if printSys
           [w,z] = damp(sys)
           zeta1
           zeta2
           omega1
           omega2
        end
    
    % Single-mode system
    else

        if ver == "rigid body"
            omega = w(3);
            zeta = z(3);
    
%             % DEBUG
%             % Check that damping ratio & nat. freq. are correct
%             omega       = sqrt(k/(m+M));
%             zeta        = b/(2*sqrt(k*(m+M)));
%             [w,z]       = damp(sys);
%             disp("Hard coded omega & zeta:")
%             omega
%             zeta
%             disp("Function calculated omega & zeta:")
%             w
%             z
    
        elseif ver == "no impedance"
            omega = sqrt(g/l);
            zeta = 0;
            omega = w(3);
            zeta = z(3);
    
        else % Fast or slow mode
            if ver == "slow"
                omega = min(w(w > 0));
                zeta = min(z(z > -1));

            elseif ver == "fast"
                omega = max(w(w > 0));
                zeta = max(z(z > -1));
            end
        end

        % If overdamped, raise flag and set other outputs to zero
        if zeta >= 1 || isnan(zeta)
            overdamped = true;
            Td = 0;
            zeta = 0;
        else
            omega_d = omega * sqrt(1 - zeta ^ 2);                           % Damped natural frequency
            Td = 2 * pi / omega_d;                                          % Damped period of vibration
            overdamped = false;
        end

        % Set multi-mode period & damping ratio outputs to zero
        Td1 = 0;
        Td2 = 0;
        zeta1 = 0;
        zeta2 = 0;

    end

end
