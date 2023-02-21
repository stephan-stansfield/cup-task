function dydt = nonlinSysEqns(t, y, u, ut, b, k)

%     global m M l g;
    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;
    
    u = interp1(ut,u,t);                                                    % Interpolate the control input (ut,u) at time t

    % Evaluate nonlinear system of equations at time t
    % y(1) = x
    % y(2) = theta
    % y(3) = x_dot
    % y(4) = theta_dot
    dydt = zeros(4, 1);
    dydt(1) = y(3);                                                         % x_dot
    dydt(2) = y(4);                                                         % theta_dot
    dydt(3) = -l * ( (g * (m + M) * tan(y(2)) + m * l * y(4) ^ 2 * ...      % x_ddot
        sin(y(2)) - k * y(1) - b * y(3) + k * u) / ...
        (m * l * cos(y(2)) ^ 2 - l * (m + M)) ) - g * tan(y(2));          
    dydt(4) = (g * (m + M) * tan(y(2)) + m * l * y(4) ^ 2 * sin(y(2)) + ... % theta_ddot  
         - k * y(1) - b * y(3) + k * u) / ...
        (m * l * cos(y(2)) - (m + M) * l / cos(y(2)));                      

end