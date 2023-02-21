function [extSys,intSys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,forwardF,intModel,impedance,printSys)

    % DEBUG
%     disp("Inside sysCreate")

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % First create system
    % Simulate system including hand impedance
    if impedance

        % External system
        A = [      0,                0,               1,      0 ;       % * [    x    ;
                   0,                0,               0,      1 ;       %      theta  ;
                -k / M,          g * m / M,        -b / M,    0 ;       %      x_dot  ;
             k / (l * M), - g / l * (m + M) / M, b / (l * M), 0 ];      %    theta_dot ]

        % Input U(t) includes feedforward force & part of impedance force
        B = [      0 ;                                                  % * U(t) = 1/k * f + b/k * xdes_dot + xdes
                   0 ;
                 k / M ;
              -k / (l * M) ];

        C = eye(4);

        D = zeros(4, 1);
        
        Arigid = [       0,       0,      1,       0 ;                  % * [    x    ;
                         0,       0,      0,       1 ;                  %      theta  ;
                    -k / (m + M), 0, -b / (m + M), 0 ;                  %      x_dot  ;    
                         0,       0,      0,       0 ];                 %    theta_dot ]
                    
        Brigid = [    0 ;                                               % * U(t)
                      0 ;
                  k / (m + M);
                      0 ];

        states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        inputs = {'u'};
        outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        extSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);
        sysRigid = ss(Arigid, Brigid, C, D, 'statename', states, ...
                'inputname', inputs, 'outputname', outputs);

        % Multi-mode, slow-mode, or fast-mode internal model
        if intModel == "full" || intModel == "slow" || intModel == "fast"
        
            intSys = extSys;
        
        % Rigid-body and no-impedance internal models
        else 

            % Simplification: Lumped mass cup and ball
            if intModel == "rigid body"
    
                A = [       0,       0,      1,       0 ;                   % * [    x    ;
                            0,       0,      0,       1 ;                   %      theta  ;
                       -k / (m + M), 0, -b / (m + M), 0 ;                   %      x_dot  ;
                            0,       0,      0,       0 ];                  %    theta_dot ]
                        
                B = [    0 ;                                                % * U(t)
                         0 ;
                      k / (m + M) ;
                         0 ];

                C = eye(4);

                D = zeros(4, 1);

            % Simplification: No hand impedance
            elseif intModel == "no impedance"

                A=[ 0,    0,   1, 0 ;                                       % * [    x    ;       
                    0,    0,   0, 1 ;                                       %      theta  ;
                    0,    0,   0, 0 ;                                       %      x_dot  ;
                    0, -g / l, 0, 0 ];                                      %    theta_dot ]

                B=[   0;                                                    % * ac (shaped x acceleration)
                      0;
                      1;
                    -1/l];

                C = eye(4);
        
                D = zeros(4, 1);
    
            end

            intSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);

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
        
        Arigid =   [0 0 1 0;
                    0 0 0 1;
                    0 0 0 0;
                    0 0 0 0];
                
        Brigid =   [0;
                    0;
                    1;
                    0];

        states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        inputs = {'u'};
        outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        extSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);
        sysRigid = ss(Arigid, Brigid, C, D, 'statename', states, ...
                'inputname', inputs, 'outputname', outputs);

        % Assign natural frequency and damping ratio to output variables
        [w,z] = damp(extSys);
        omega = w(3);
        zeta = z(3);
        Td = 2*pi/omega;
        overdamped = false;
        
    end
    
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
        pzmap(extSys)

    end

    %%%%%%%
    % Get natural frequency and damping ratio of internal model
    [w,z] = damp(intSys);

    % Compute internal model natural frequencies and damping ratios
    if intModel == "full"

        % Assign natural frequency and damping ratio of two modes to variables
        w1 = w(1);
        w2 = w(3);
        z1 = z(1);
        z2 = z(3);

        % Input shaper using LOWER frequency
        zeta1 = z1;
        omega1 = w1;

        % Input shaper using HIGHER frequency
        zeta2 = z2;
        omega2 = w2;

        % Check if either mode is overdamped
        if zeta1 >= 1 || zeta2 >= 1

            % Raise flag and set other outputs to zero
            overdamped = true;
            Td1 = 0;
            Td2 = 0;
            zeta1 = 0;
            zeta2 = 0;

        else

            omega_d1 = omega1 * sqrt(1 - zeta1 ^ 2);                        % Damped natural frequency, lower mode (rad/s)
            Td1 = 2 * pi / omega_d1;                                        % Damped period of vibration, lower mode (s)

            omega_d2 = omega2 * sqrt(1 - zeta2 ^ 2);                        % Damped natural frequency, higher mode (rad/s)
            Td2 = 2 * pi / omega_d2;                                        % Damped period of vibration, higher mode (s)
            overdamped = false;

        end
        
        % Set single-mode period & damping ratio outputs to zero
        Td = 0;
        zeta = 0;
        
        if printSys
           [w,z] = damp(extSys)
           zeta1
           zeta2
           omega1
           omega2
        end
    
    else

        % Single-mode internal models
        if intModel == "rigid body" || intModel == "no impedance"

            omega = w(3);
            zeta = z(3);
    
%             % DEBUG
%             % Check that damping ratio & nat. freq. are correct
%             % For rigid body:
%             omega       = sqrt(k/(m+M));
%             zeta        = b/(2*sqrt(k*(m+M)));
%             [w,z]       = damp(intSys);
%             disp("Hard coded omega & zeta:")
%             omega
%             zeta
%             disp("Function calculated omega & zeta:")
%             w
%             z

%             % For no impedance:
%             omega = sqrt(g/l);
%             zeta = 0;
    
        else 
            
            % Fast or slow mode
            if intModel == "slow"

                [omega, min_ind] = min(w(w > 0));
                zeta_temp = z(w > 0);
                zeta = zeta_temp(min_ind);
                m_mod = M;

%                 % DEBUG
%                 disp('In sysCreate')
%                 w
%                 min_ind
%                 z
%                 omega
%                 zeta

            elseif intModel == "fast"

                [omega, max_ind] = max(w(w > 0));
                zeta_temp = z(w > 0);
                zeta = zeta_temp(max_ind);
                m_mod = m;

%                 % DEBUG
%                 disp('In sysCreate')
%                 w
%                 max_ind
%                 z
%                 omega
%                 zeta

            end
            
            % Calculate effective modal stiffness and damping using larger
            % mass (cart) for slow mode and smaller mass (ball) for fast mode
            k_mod = m_mod * omega ^ 2;
            b_mod = 2 * m_mod * zeta * omega;

            % Create new 2nd-order system using effective modal parameters
            A = [       0,       0,        1,       0 ;                     % * [    x 
                        0,       0,        0,       1 ;                     %      theta
                 -k_mod / m_mod, 0, -b_mod / m_mod, 0 ;                     %      x_dot    
                        0,       0,        0,       0 ];                    %    theta_dot ]
                        
            B = [      0 ;                                                  % * U(t)
                       0 ;
                 k_mod / m_mod ;
                       0 ];

            C = eye(4);

            D = zeros(4, 1);

            states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
            inputs = {'u'};
            outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
            intSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, ...
                'outputname', outputs);

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
