function [sys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = ...
        sysCreate(b,k,forwardF,ver,impedance,printSys)

    % DEBUG
%     disp("Inside sysCreate")

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % Simulate system including hand impedance. Linearized 1D pendulum
    % below cart w/ stiffness & damping (used in Maurice et al, 2018)
    if impedance
        
        % Controller includes feedforward force input term
        if forwardF
            
            A = [      0,                0,               1,      0 ;       % * [    x
                       0,                0,               0,      1 ;       %      theta
                    -k / M,          g * m / M,        -b / M,    0 ;       %      x_dot
                 k / (l * M), - g / l * (m + M) / M, b / (l * M), 0 ];      %    theta_dot ]

            % Input U(t) includes feedforward force and half of impedance
            % force
            B = [      0 ;                                                  % * U(t)
                       0 ;
                     k / M ;
                  -k / (l * M) ];

            C = eye(4);

            D = zeros(4, 1);

            states = {'x', 'theta', 'x_{dot}', 'theta_{dot}'};
            
            Arigid = [       0,       0,      1,       0 ;                  % * [    x 
                             0,       0,      0,       1 ;                  %      theta
                        -k / (m + M), 0, -b / (m + M), 0 ;                  %      x_dot    
                             0,       0,      0,       0 ];                 %    theta_dot ]
                        
           Brigid = [    0 ;                                                % * U(t)
                         0 ;
                     k / (m + M);
                         0 ];

        % No feedforward force input term; forwardF = false
        else
%             %{
%             A=[    0       0            0          0        0 ;             % *[x_des, x, theta, x_dot, theta_dot]'
%                    0       0            0          1        0 ;        
%                    0       0            0          0        1 ;
%                   k/M    -k/M         m*g/M      -b/M       0 ;
%                -k/(l*M) k/(l*M) -g*(m+M)/(l*M)  b/(l*M)     0 ];
% 
%             % Input "u" is desired velocity
%             B=[   1    ;                                                    % *x_des_dot
%                   0    ;                 
%                   0    ;
%                  b/M   ;
%                -b/(l*M)];
%            
%            % Note: added 5th output, xdes, for pendulum lock code.
%            C=[ 0      1     0      0      0;                                % *[x_des, x, theta, x_dot, theta_dot]'
%                0      0 360/(2*pi) 0      0;
%                0      0      0     1      0;
%                0      0      0     0   360/(2*pi);
%                1      0      0     0      0];
% 
%             D=[0;                                                           % *x_des_dot
%                0;
%                0;
%                0;
%                0];
% 
%             states={'x_{des}' 'x' 'theta' 'x_{dot}' 'theta_{dot}'};
%             
%             Arigid = [     0       0     0    0      0;                     % *[x_des, x, theta, x_dot, theta_dot]'
%                            0       0     0    1      0;
%                            0       0     0    0      1;
%                         k/(m+M) -k/(m+M) 0 -b/(m+M)  0;
%                            0       0     0    0      0];
%                     
%            Brigid = [   1;                                                  % *x_des_dot
%                         0;
%                         0;
%                       b/(m+M);
%                         0];
%             %}
% 
%             % TESTING %
%             % Try with x''des input instead. Does this change anything?
%             A=  [    0       0            0         1         0    0 ;      % *[x_des, x, theta, x_des_dot, x_dot, theta_dot]'
%                      0       0            0         0         1    0 ;
%                      0       0            0         0         0    1 ;
%                      0       0            0         0         0    0 ;
%                     k/M    -k/M         m*g/M      b/M      -b/M   0 ;
%                  -k/(l*M) k/(l*M) -g*(m+M)/(l*M) -b/(l*M)  b/(l*M) 0 ];
% 
%             % Input "u" is desired ACCELERATION
%             B = [   0    ;                                                  % *x_des_dotdot
%                     0    ;                 
%                     0    ;
%                     1    ;
%                     0   ;
%                     0];
%            
%            % Note: added 5th output, xdes, for pendulum lock code.
%            C=[ 0      1     0      0    0      0;                           % *[x_des, x, theta, x_des_dot, x_dot, theta_dot]'
%                0      0 360/(2*pi) 0    0      0;
%                0      0      0     0    1      0;
%                0      0      0     0    0   360/(2*pi);
%                1      0      0     0    0      0];
% 
%             D=[0;                                                           % *x_des_dotdot
%                0;
%                0;
%                0;
%                0];
% 
%             states={'x_{des}' 'x' 'theta' 'x_{dot}_{des}' 'x_{dot}' 'theta_{dot}'};
%             
%             Arigid = [     0       0     0      1      0      0;                     % *[x_des, x, theta, x_des_dot, x_dot, theta_dot]'
%                            0       0     0      0      1      0;
%                            0       0     0      0      0      1;
%                            0       0     0      0      0      0;
%                         k/(m+M) -k/(m+M) 0  b/(m+M) -b/(m+M)  0;
%                            0       0     0      0      0      0];
%                     
%            Brigid = [   0;                                                  % *x_des_dotdot
%                         0;
%                         0;
%                         1;
%                         0;
%                         0];
        end

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

        % Input shaping based on simplified internal model
        if ver ~= ""

            % Simplification: Lumped mass cup and ball
            if ver == "rigid body"
                omega       = sqrt(k/(m+M));
                zeta        = b/(2*sqrt(k*(m+M)));

            % Simplification: No hand impedance
            elseif ver == "no impedance"
                omega       = sqrt(g/l);
                zeta        = 0;

            else
                [w,z]       = damp(sys);

                % Simplification: Slower mode of full system
                if ver == "slow"
                    omega   = min(w(w>0));
                    zeta    = min(z(z>-1));

                % Simplification: Faster mode of full system
                elseif ver == "fast"
                    omega   = max(w(w>0));
                    zeta    = max(z(z>-1));
                end
            end

            if zeta >= 1 || isnan(zeta)
                % If overdamped, raise flag and set other outputs to zero
                overdamped  = true;
                Td          = 0;
                zeta        = 0;
            else
                omega_d     = omega*sqrt(1-zeta^2);                         % Damped natural frequency
                Td          = 2*pi/omega_d;                                 % Damped period of vibration
                overdamped  = false;
            end
            
            
            % Assign null values to alternate version outputs
            Td1 = 0;
            Td2 = 0;
            zeta1 = 0;
            zeta2 = 0;

        % Multi-mode input shaping based on full internal model
        else
            % Assign natural frequency and damping ratio of two modes to variables
            [w,z] = damp(sys);
%             if forwardF
                w1 = w(1);
                w2 = w(3);
                z1 = z(1);
                z2 = z(3);

%             else
%                 w1 = w(2);
%                 w2 = w(4);
%                 z1 = z(2);
%                 z2 = z(4);
                
%                 % TESTING %
%                 w1 = w(3);
%                 w2 = w(5);
%                 z1 = z(3);
%                 z2 = z(5);
%             end

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
                omega_d1 = omega1*sqrt(1-zeta1^2);                          % Damped natural frequency, lower mode (rad/s)
                Td1 = 2*pi/omega_d1;                                        % Damped period of vibration, lower mode (s)
                overdamped = false;
            end

            % Input shaper using HIGHER frequency
            zeta2 = z2;
            omega2 = w2;

            % If first mode isn't overdamped, check second mode
            if not(overdamped)
                if zeta2>=1
                    % If overdamped, raise flag and set other outputs to zero
                    overdamped = true;
                    Td1 = 0;
                    Td2 = 0;
                    zeta1 = 0;
                    zeta2 = 0;
                else
                    omega_d2=omega2*sqrt(1-zeta2^2);                        % Damped natural frequency, higher mode (rad/s)
                    Td2 = 2*pi/omega_d2;                                    % Damped period of vibration, higher mode (s)
                    overdamped = false;
                end
            end
            
            % Assign null values to alternate version outputs
            Td = 0;
            zeta = 0;
            
            if printSys
               [w,z] = damp(sys)
               zeta1
               zeta2
               omega1
               omega2
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

        inputs={'u'};

        outputs={'x' 'theta' 'x_{dot}' 'theta_{dot}'};

        sys=ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
        sysRigid=ss(Arigid,Brigid,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

        % Assign natural frequency and damping ratio to output variables
        [w,z] = damp(sys);
        omega = w(3);
        zeta = z(3);
        Td = 2*pi/omega;
        overdamped = false;
        
        % Assign null values to alternate version outputs
        Td1 = 0;
        Td2 = 0;
        zeta1 = 0;
        zeta2 = 0;
        
    end

end
