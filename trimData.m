function [pos,theta,vel,omega,acc,alpha,tdes,tc,tRaw,velRaw,accRaw] = ...
        trimData(pos,theta,vel,omega,acc,alpha,t,st,start,stop)
% TRIMDATA
%
% Trims cup task experimental data from HapticMaster given particular cart
% position and velocity characteristics.

    % Convert ball units to degrees
    theta = rad2deg(theta);
    omega = rad2deg(omega);
    alpha = rad2deg(alpha);
    
    % Offset cart position so it starts at x = 0
    pos = pos + 0.125;
    
    % Trim vectors to specified start and stop indices
    if stop < length(t)
        t       = t(start:stop);
        pos     = pos(start:stop);
        theta   = theta(start:stop);
        vel     = vel(start:stop);
        omega   = omega(start:stop);
        acc     = acc(start:stop);
        alpha   = alpha(start:stop);
    else
        pos     = pos(start:end);
        theta   = theta(start:end);
        vel     = vel(start:end);
        omega   = omega(start:end);
        acc     = acc(start:end);
        alpha   = alpha(start:end);
        
        pos(length(t)) = pos(end);
        vel(length(t)) = vel(end);
        acc(length(t)) = acc(end);
        theta(length(t)) = theta(end);
        omega(length(t)) = omega(end);
        alpha(length(t)) = alpha(end);
    end

    % Shift experimental time vector and create linear time vector
    t       = t - t(1);
    tdes    = t(end) - t(1);
    tc      = 0:st:tdes;
    
    % Return uninterpolated values of interest
    tRaw   = t;
    velRaw = vel;
    accRaw = acc;
    
    % Interpolate data to match simulated sample points
    pos     = interp1(t,pos,tc);
    theta   = interp1(t,theta,tc);
    vel     = interp1(t,vel,tc);
    omega   = interp1(t,omega,tc);
    acc     = interp1(t,acc,tc);
    alpha   = interp1(t,alpha,tc);
    
end