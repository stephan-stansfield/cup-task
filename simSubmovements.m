function [output, vc] = simSubmovements(sys,sysRigid,tExp,tDesSim,st,xEnd,vStart,...
        D1,tf1,ti2,pendIndex,simVersion,forwardF,impedance)
%% Generate desired input using minimum-jerk profile (min jerk of hand)

% Reset flag for invalid trials
flag = false;

% Run first-pass simulation
[output, startIndex] = sim();

if flag
    % If first run produces a rejection, take the output and return to
    % the optimization script without running a second time
    return
end

% Add amount of time that was trimmed to pendulum release time
pendIndex = pendIndex + startIndex;

% Simulate again with lengthened duration
[output, ~] = sim();

function [output, startIndex] = sim()

    t = 0:st:tDesSim;
    lenT = length(t);

    xi1 = 0;                                                                % Initial position for 1st submovement
    xf1 = D1;                                                               % Final position for 1st submovement

    xi2 = D1;                                                               % Initial position for 2nd submovement
    xf2 = xEnd;                                                             % Final position for 2nd submovement

    tf2 = tDesSim;                                                          % Final time for 2nd submovement

    % Generate 1st minimum-jerk trajectory submovement
    t1 = 0:st:tf1;
    a1 = (xf1-xi1)*(60*(t1/tf1) - 180*(t1/tf1).^2 ...                       % Acceleration profile for 1st submovement
        + 120*(t1/tf1).^3)/(tf1^2);   
    v1 = (xf1-xi1)*(30*(t1/tf1).^2 - 60*(t1/tf1).^3 + 30*(t1/tf1).^4)/tf1;

    % Pad end of submovement 1 with zeros
    padlength   = length(t) - length(a1);
    a1          = [a1, zeros(1,padlength)];
    v1          = [v1, zeros(1,padlength)];

    % Generate 2nd minimum-jerk trajectory submovement
    t2 = ti2:st:tf2;
    a2 = (xf2-xi2)*(60*((t2-ti2)/(tf2-ti2)) - 180*((t2-ti2)/(tf2-ti2)).^2 ...   % Acceleration profile for 2nd submovement
        + 120*((t2-ti2)/(tf2-ti2)).^3)/((tf2-ti2)^2);      
    v2 = (xf2-xi2)*(30*((t2-ti2)/(tf2-ti2)).^2 - 60*((t2-ti2)/(tf2-ti2)).^3 ...
        + 30*((t2-ti2)/(tf2-ti2)).^4)/(tf2-ti2);

    % Pad beginning of submovement 2 with zeros
    padlength   = length(t) - length(a2);
    a2          = [zeros(1,padlength), a2];
    v2          = [zeros(1,padlength), v2];

    % Add submovements
    atot        = a1 + a2;
    vtot        = v1 + v2;
    vc          = [v1; v2];                                                 % Stack submovements in array for plotting
    
    if simVersion == "linear"
        % Simulate system response using acceleration profile as simulation input.
        % Store results in output variable.
        output      = lsim(sys,atot,t);

        % Assign outputs to variables
        pos_sim     = output(:,1);
        theta_sim   = output(:,2);
        vel_sim     = output(:,3);
        omega_sim   = output(:,4);

    elseif simVersion == "pendLock"
        % Reject profiles that end before the pendulum is released
        if pendIndex >= lenT
            startIndex = 0;
            output = 100000*ones(lenT,6);
            flag = true;
            return
        end   

        % Divide time and input vectors into two sections: before and
        % after time that pendulum starts to move (when cart velocity
        % exceeds 0.1 m/s)
        t1 = t(1:pendIndex);
        t2 = t(pendIndex+1:end)-t(pendIndex+1);                             % Shift value so first t = 0
        aPrePend = atot(1:pendIndex);
        aPostPend = atot(pendIndex+1:end);

        % Simulate rigid-body system for first portion of movement
        output1  = lsim(sysRigid,aPrePend,t1);
        
        % Capture kinematic values at end of first section, X0
        if impedance
            if forwardF
                if length(output1(1,:)) ~= 4
                    X0 = [output1(end,5), output1(end,1),...                % xdes, x,
                        output1(end,2)*2*pi/360, vtot(pendIndex), ...       % theta (rad), xdes_dot = vtot,
                        output1(end,3), output1(end,4)*2*pi/360];           % x_dot, theta_dot (rad/s)
                else
                    X0 = [output1(end,1);                                   % x
                          output1(end,2)*2*pi/360;                          % theta (rad)
                          output1(end,3);                                   % x_dot
                          output1(end,4)*2*pi/360];                         % theta_dot (rad/s)
                end
            else
                X0 = [output1(end,5), output1(end,1),...                    % xdes, x,
                    output1(end,2)*2*pi/360, vtot(pendIndex),...            % theta (rad), xdes_dot = vc,
                    output1(end,3), output1(end,4)*2*pi/360];               % x_dot, theta_dot (rad/s)
            end
        else
            X0 = [output1(end,1), output1(end,2)*2*pi/360,...               % x, theta (rad), 
                output1(end,3), output1(end,4)*2*pi/360];                   % x_dot, theta_dot (rad/s)
        end
        
        % When pendulum is released, simulate full system for the
        % remainder of the move duration. Use values at end of previous
        % simulation, X0, as initial conditions for this simulation
        output2 = lsim(sys,aPostPend,t2,X0);
        
        % Concatenate outputs
        output = [output1; output2];

        % Assign outputs to variables
        pos_sim     = output(:,1);
        theta_sim   = output(:,2);
        vel_sim     = output(:,3);
        omega_sim   = output(:,4);
        
    end

    % Take approximative derivative of cart velocity to get cart
    % acceleration profile. Take central difference at all interior points,
    % forward difference at first point, and backward difference at last
    % point.
    acc_sim         = zeros(length(vel_sim),1);                             % Initialize cart acceleration vector
    acc_sim(1)      = (vel_sim(2) - vel_sim(1))/st;                         % Forward difference
    acc_sim(end)    = (vel_sim(end) - vel_sim(end-1))/st;                   % Backward difference
    for i = 2:1:length(acc_sim)-1
        acc_sim(i)  = (vel_sim(i+1) - vel_sim(i-1))/(2*st);                 % Central difference
    end

    % Take approximate derivative of ball velocity in same way as cart
    % to get ball angular acceleration profile.
    alpha_sim           = zeros(length(omega_sim),1);                       % Initialize ball acceleration vector
    alpha_sim(1)        = (omega_sim(2) - omega_sim(1))/st;                 % Forward difference
    alpha_sim(end)      = (omega_sim(end) - omega_sim(end-1))/st;           % Backward difference
    for i = 2:1:length(omega_sim)-1
        alpha_sim(i)    = (omega_sim(i+1) - omega_sim(i-1))/(2*st);         % Central difference
    end
        
    % Trim simulated output at first point where cart velocity exceeds 
    % experimental start velocity
    startIndex      = find(vel_sim > vStart,1);
    if ~isempty(startIndex)
        pos_sim     = pos_sim(startIndex:end);
        theta_sim   = theta_sim(startIndex:end);
        vel_sim     = vel_sim(startIndex:end);
        omega_sim   = omega_sim(startIndex:end);
        acc_sim     = acc_sim(startIndex:end);
        alpha_sim   = alpha_sim(startIndex:end);
        vc          = vc(:,startIndex:end);
    else
        startIndex  = 0;
    end

    % If simulated profile is longer than experimental profile, trim to
    % same length
    lenExp = length(0:st:tExp);
    if length(vel_sim) > lenExp
        pos_sim     = pos_sim(1:lenExp);
        theta_sim   = theta_sim(1:lenExp);
        vel_sim     = vel_sim(1:lenExp);
        omega_sim   = omega_sim(1:lenExp);
        acc_sim     = acc_sim(1:lenExp);
        alpha_sim   = alpha_sim(1:lenExp);
        vc          = vc(:,1:lenExp);
    end

    % Re-assign trimmed kinematic profiles to output array
    output = [pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim];

end

end
