function [extSys,intSys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,intModel,impedance,printSys)
% SYSCREATE
% Creates a linear state-space system of the form:
%   dx = A*x + B*u
%   y = C*x + D*u
% representing some model of the cart-and-pendulum system used in the cup
% task. This function is used to represent both the "full" system of the
% cart-and-pendulum coupled to a human's upper limb's impedance, as well as
% simplified representations of that system that humans may use as an
% "internal model".

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % Simulate system including hand impedance
    if impedance
        A = [      0,                0,               1,      0 ;           % * [    x    ;
                   0,                0,               0,      1 ;           %      theta  ;
                -k / M,          g * m / M,        -b / M,    0 ;           %      x_dot  ;
             k / (l * M), - g / l * (m + M) / M, b / (l * M), 0 ];          %    theta_dot ]

        % Input U(t) includes feedforward force & part of impedance force
        B = [      0 ;                                                      % * U(t) = f + b*xdes_dot + k*xdes
                   0 ;
                 1 / M ;
              -1 / (l * M) ];

        C = eye(4);

        D = zeros(4, 1);
        
        % "rigid" state and input matrices are used to recreate
        % experimental bug where pendulum was locked at start of trials
        Arigid = [       0,       0,      1,       0 ;                      % * [    x    ;
                         0,       0,      0,       1 ;                      %      theta  ;
                    -k / (m + M), 0, -b / (m + M), 0 ;                      %      x_dot  ;    
                         0,       0,      0,       0 ];                     %    theta_dot ]
                    
        Brigid = [    0 ;                                                   % * U(t) = f + b*xdes_dot + k*xdes
                      0 ;
                  1 / (m + M);
                      0 ];

        % Create "external" system including cart, pendulum, and coupled
        % hand impedance
        states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        inputs = {'u'};
        outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
        extSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs,...
            'outputname', outputs);
        sysRigid = ss(Arigid, Brigid, C, D, 'statename', states, ...
                'inputname', inputs, 'outputname', outputs);

        if intModel == "full" || intModel == "slow" || intModel == "fast"
            intSys = extSys;
        else 
            if intModel == "rigid body"
    
                A = [       0,       0,      1,       0 ;                   % * [    x    ;
                            0,       0,      0,       1 ;                   %      theta  ;
                       -k / (m + M), 0, -b / (m + M), 0 ;                   %      x_dot  ;
                            0,       0,      0,       0 ];                  %    theta_dot ]
                        
                B = [    0 ;                                                % * U(t) = f + b*xdes_dot + k*xdes
                         0 ;
                      1 / (m + M) ;
                         0 ];

                C = eye(4);

                D = zeros(4, 1);

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

            % Create "internal" system based on simplified model
            intSys = ss(A,B,C,D, 'statename', states, 'inputname', inputs,...
                'outputname', outputs);

        end
    
    else
        % Simulate system without hand impedance
        A=[0   0  1 0;                                                      % *[x theta x_dot theta_dot]'
           0   0  0 1;
           0   0  0 0;
           0 -g/l 0 0];

        % Single input "u" is desired acceleration
        B=[   0;                                                            % *xdes_ddot'
              0;
              1;
            -1/l];

        % 4 outputs: cart pos, ball angle, cart vel, & ball angular vel
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
        intSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);
        extSys = ss(A, B, C, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);
        sysRigid = ss(Arigid, Brigid, C, D, 'statename', states, ...
                'inputname', inputs, 'outputname', outputs);

        % Assign natural frequency and damping ratio to output variables
        [w,z] = damp(extSys);
        omega = w(3);
        zeta = z(3);
        Td = 2 * pi / omega;
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
    % Compute internal model's natural frequencies and damping ratios
    [w,z] = damp(intSys);

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
            % Lower mode damped natural frequency (rad/s) and damped period
            % of vibration (s)
            omega_d1 = omega1 * sqrt(1 - zeta1 ^ 2);
            Td1 = 2 * pi / omega_d1;

            % Higher mode damped natural frequency (rad/s) and damped
            % period of vibration (s)
            omega_d2 = omega2 * sqrt(1 - zeta2 ^ 2);
            Td2 = 2 * pi / omega_d2;
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
    
        else 
            
            % Fast or slow mode
            if intModel == "slow"
                [omega, min_ind] = min(w(w > 0));
                zeta_temp = z(w > 0);
                zeta = zeta_temp(min_ind);
                m_mod = M;
            elseif intModel == "fast"
                [omega, max_ind] = max(w(w > 0));
                zeta_temp = z(w > 0);
                zeta = zeta_temp(max_ind);
                m_mod = m;
            end
            
            % Calculate effective modal stiffness and damping using larger
            % mass (cart) for slow mode and smaller mass (ball) for fast
            % mode
            k_mod = m_mod * omega ^ 2;
            b_mod = 2 * m_mod * zeta * omega;

            % Create new 2nd-order system using effective modal parameters
            A = [       0,       0,        1,       0 ;                     % * [    x 
                        0,       0,        0,       1 ;                     %      theta
                 -k_mod / m_mod, 0, -b_mod / m_mod, 0 ;                     %      x_dot    
                        0,       0,        0,       0 ];                    %    theta_dot ]
                        
            B = [      0 ;                                                  % * U(t)
                       0 ;
                 1 / m_mod ;
                       0 ];

            C = eye(4);

            D = zeros(4, 1);

            states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
            inputs = {'u'};
            outputs = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
            intSys = ss(A,B,C,D, 'statename', states, 'inputname', inputs,...
                'outputname', outputs);

        end

        % If overdamped, raise flag and set other outputs to zero
        if zeta >= 1 || isnan(zeta)
            overdamped = true;
            Td = 0;
            zeta = 0;
        else
            % Damped natural frequency and damped period of vibration
            omega_d = omega * sqrt(1 - zeta ^ 2);
            Td = 2 * pi / omega_d;
            overdamped = false;
        end

        % Set multi-mode period & damping ratio outputs to zero
        Td1 = 0;
        Td2 = 0;
        zeta1 = 0;
        zeta2 = 0;

    end
end