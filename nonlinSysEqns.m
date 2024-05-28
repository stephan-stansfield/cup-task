function dydt = nonlinSysEqns(t, y, u, ut, b, k)
% NONLINSYSEQNS
% Evaluates nonlinear system of equations of cart-and-pendulum object
% coupled to stiffness and damping of upper limb at time t given control
% input u.

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;
    
    % Interpolate the control input (ut,u) at time t
    u = interp1(ut,u,t);

    % Evaluate nonlinear system of equations at time t
    % y(1) = x
    % y(2) = theta
    % y(3) = x_dot
    % y(4) = theta_dot
    dydt = zeros(4, 1);
    dydt(1) = y(3);
    dydt(2) = y(4);
    dydt(3) = -l * ( (g * (m + M) * tan(y(2)) + m * l * y(4) ^ 2 * ...      
        sin(y(2)) - k * y(1) - b * y(3) + u) / ...
        (m * l * cos(y(2)) ^ 2 - l * (m + M)) ) - g * tan(y(2));          
    dydt(4) = (g * (m + M) * tan(y(2)) + m * l * y(4) ^ 2 * sin(y(2)) ...  
         - k * y(1) - b * y(3) + u) / ...
        (m * l * cos(y(2)) - (m + M) * l / cos(y(2)));                      

end