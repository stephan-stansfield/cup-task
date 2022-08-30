function [sys,Td,zeta] = corrSysCreate(b,k)

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % Simulate system including hand impedance. Linearized 1D pendulum
    % below cart w/ stiffness & damping (used in Maurice et al, 2018)
            
    A=[    0       0             0            1        0   0 ;              % * [x_des  
           0       0             0            0        1   0 ;              %    x
           0       0             0            0        0   1 ;              %    theta
           0       0             0            0        0   0 ;              %    x_des_dot
          k/M    -k/M          m*g/M         b/M    -b/M   0 ;              %    x_dot
       -k/(l*M) k/(l*M) -(g*(m+M))/(l*M)  -b/(l*M) b/(l*M) 0 ];             %    theta_dot]

    % Single input "u" is desired acceleration
    B=[      0      ;                                                       % * x_des_dotdot
             0      ;
             0      ;
             1      ;
          (m+M)/M   ;
       -(m+M)/(l*M) ];

    % 5 outputs: cart position, ball angle, cart vel, & ball
    % angular velocity. Ball angle converted from radians to deg.
    % 5th output is xdes for pendulum lock code to fully "hand off"
    % the state between simulations when pendulum is released.
    C=[0      1     0      0    0     0;                                    % * [x_des  
       0      0 360/(2*pi) 0    0     0;                                    %    x
       0      0     0      0    1     0;                                    %    theta
       0      0     0      0    0 360/(2*pi);                               %    x_des_dot
       1      0     0      0    0     0];                                   %    x_dot
                                                                            %    theta_dot]

    D=[0;                                                                   % * x_des_dotdot
       0;
       0;
       0;
       0];

    states={'x_{des}' 'x' 'theta' 'x_{dot}_{des}' 'x_{dot}' 'theta_{dot}'};

%     Arigid = [     0       0     0    1        0    0;                      % *[x_des, x, theta, x_des_dot, x_dot, theta_dot]'
%                    0       0     0    0        1    0;
%                    0       0     0    0        0    1;
%                    0       0     0    0        0    0;
%                 k/(m+M) -k/(m+M) 0 b/(m+M) -b/(m+M) 0;
%                    0       0     0    0        0    0];
% 
%    Brigid = [0;                                                             % *x_des_dotdot
%              0;
%              0;
%              1;
%              1;
%              0];

    inputs={'u'};
    outputs={'x' 'theta' 'x_{dot}' 'theta_{dot}' 'x_{des}'};
    sys=ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
%     sysRigid=ss(Arigid,Brigid,C,D,'statename',states,'inputname',inputs,'outputname',outputs);

    % Simplification: Lumped mass cup and ball
    omega       = sqrt(k/(m+M));
    zeta        = b/(2*sqrt(k*(m+M)));

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

end
