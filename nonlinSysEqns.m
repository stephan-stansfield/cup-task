function dydt = nonlinSysEqns(t,y,u,ut,b,k)

%     global m M l g;
    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;
    
    u = interp1(ut,u,t);                                                    % Interpolate the control input (ut,u) at time t

    % Evaluate nonlinear system of equations at time t
    dydt = zeros(6,1);
    dydt(1) = y(4);                                                         % x_des_dot
    dydt(2) = y(5);                                                         % x_dot
    dydt(3) = y(6);                                                         % th_dot
    dydt(4) = u;                                                            % x_des_ddot
    dydt(5) = -l*((g*(m+M)*tan(y(3)) + m*l*y(6)^2*sin(y(3))+(M+m)*u-...     % x_ddot
        k*(y(2)-y(1))-b*(y(5)-y(4)))/(m*l*(cos(y(3)))^2-l*(m+M)))-g*tan(y(3));    
    dydt(6) = (g*(m+M)*tan(y(3)) + m*l*y(3)^2*sin(y(3)) + (M+m)*u-...       % th_ddot
        k*(y(2)-y(1))-b*(y(5)-y(4)))/(m*l*cos(y(3)) - (m+M)*l/cos(y(3)));   

end