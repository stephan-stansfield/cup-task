function [output, dur_corr, vc] = simInputShape(b,k,intModel,extSys,...
        sysRigid,Td1,Td2,zeta1,zeta2,tExp,tDesSim,xEnd,vStart,st,forwardF,...
        simVersion,modes,pendIndex,impedance)
% SIMINPUTSHAPE
%
% This function generates 2 or 4 input shaping impulses based on a system
% with 1 or 2 modes of oscillation, respectively. The impulses will be
% symmetric if the system is lossless (zeta = 0). If the system is
% underdamped (0 < zeta < 1), the impulses will be asymmetrical. If the
% system is critically or overdamped (zeta >= 1), input shaping will not
% work.

    G = globalData();
    m = G.m;
    M = G.M;
    g = G.g;
    
    % Initialize flag for invalid trials
    invalidFlag = false;
    
    % Run first-pass simulation to determine pendulum release time
    [output, startIndex, dur_corr] = sim();
    
    % If first run produces a rejection, take the output and return to
    % the optimization function without running a second time
    if invalidFlag
        return
    end
             
    % Add duration that was trimmed from start to pendulum release time
    pendIndex = pendIndex + startIndex;
    
    % Simulate again with lengthened duration
    [output, ~, dur_corr] = sim();
    
    function [output, startIndex, dur_corr] = sim()

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (i): Generate Input Shaping impulses
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % 1 mode (2 input shaping impulses)
        if modes == 1

            % Generate impulse times
            % @tf: duration of pre-convolved profile
            % @t: time vector pre-convolution
            % @tc: time vector post-convolution
            t1 = 0;                                                            
            t2 = Td1 / 2;
            t2 = round(t2, mydigits(st));
            tf = tDesSim - t2;
            t  = 0:st:tf;
            tc = 0:st:tDesSim;
            lentc = length(tc);

            if tf <= 0
                startIndex = 0;
                output = 100000*ones(length(tc), 6);
                vc = [];
                dur_corr = 0;
                invalidFlag = true;
                return
            end

            % Generate impulse amplitudes
            K = exp(-zeta1 * pi / sqrt(1 - zeta1 ^ 2));     
            A1 = 1 / (1 + K);
            A2 = K / (1 + K);
            
            % Create vector of impulses. Assign amplitudes A1 and A2 to
            % times t1 and t2, respectively. Trim trailing zeros.
            pulses = zeros(1, length(t));                           
            pulses(int16(1 + t1 / st)) = A1;
            pulses(int16(1 + t2 / st)) = A2;
            pulsetimes = find(pulses);
            pulses = pulses(1:pulsetimes(end));
        
        % 2 modes (4 input shaping impulses)
        elseif modes == 2 
            % Generate impulse times
            t1 = 0;
            t21 = Td1 / 2;
            t22 = Td2 / 2;
            toffset = t21 + t22;
            tf = tDesSim - toffset;
            t = 0:st:tf;
            tc = 0:st:tDesSim;
            lentc = length(tc);

            if tf <= 0
                startIndex = 0;
                output = 100000*ones(length(tc), 6);
                vc = [];
                dur_corr = 0;
                invalidFlag = true;
                return
            end

            % Generate first mode impulse amplitudes
            K1 = exp(-zeta1 * pi / sqrt(1 - zeta1 ^ 2));     
            A11 = 1 / (1 + K1);
            A21 = K1 / (1 + K1);

            % Generate second mode impulse amplitudes
            K2 = exp(-zeta2 * pi / sqrt(1 - zeta2 ^ 2));     
            A12 = 1 / (1 + K2);                        
            A22 = K2 / (1 + K2);             

            % Create vector of first mode impulses
            t21 = round(t21, mydigits(st));
            imp1 = zeros(1, length(t));
            imp1(int16(1 + t1 / st)) = A11;
            imp1(int16(1 + t21 / st)) = A21;

            % Create vector of second mode impulses
            t22 = round(t22, mydigits(st));
            imp2 = zeros(1, length(t));
            imp2(int16(1 + t1 / st)) = A12;
            imp2(int16(1 + t22 / st)) = A22;

            % Convolve the two sets of impulses and trim trailing zeros.
            pulses = conv(imp1, imp2);
            pulsetimes = find(pulses);
            pulses = pulses(1:pulsetimes(end));

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (ii): Shape desired input
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Minimum-jerk profiles
        v = xEnd * (30*(t/tf).^2 - 60*(t/tf).^3 + 30*(t/tf).^4) / tf;
        a = xEnd * (60*(t/tf) - 180*(t/tf).^2 + 120*(t/tf).^3) / (tf ^ 2);

        % Convolve min jerk profiles with input shaping impulses
        vc = conv(v, pulses);
        ac = conv(a, pulses);

        % Make time and kinematics arrays same length if rounding error
        if length(vc) > lentc % Trim vectors after desired time  
            vc = vc(1:lentc);
            ac = ac(1:lentc);
        else % Pad end of array with zeros
            vc(lentc) = 0;
            ac(lentc) = 0;
        end

        % Integrate to get convolved position profile. Correct for cumtrapz
        % calculating 0 for last index
        xc = cumtrapz(tc, vc);
        xc(end) = xc(end - 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iii): Generate control input
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if forwardF
            % If using feedforward force, simulate the response of the
            % internal model to the shaped input and use this as the
            % "desired" trajectory in the inverse dynamics

            % Create dynamic system of internal model
            [~,intSys,~,~,~,~,~,~,~,~] = sysCreate(b, k, intModel,...
                impedance, false);

            % Choose input by internal model
            if intModel == "full" || intModel == "rigid body"
                u = b * vc + k * xc ;
            elseif intModel == "no impedance"
                u = ac;
            elseif intModel == "slow" || intModel == "fast"
                [w, z] = damp(intSys);
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

                k_mod = m_mod * omega ^ 2;
                b_mod = 2 * m_mod * zeta * omega;
                u = b_mod * vc + k_mod * xc ;
            end

            % Simulate response of internal model when given shaped input
            int_mdl_output = lsim(intSys, u, tc);

            % Get desired position and velocity profiles (output of
            % internal model simulation)
            x_des = int_mdl_output(:, 1);
            v_des = int_mdl_output(:, 3);

            % Differentiate desired velocity, taking central difference at
            % interior points and forward and backward differences at first
            % and last points to get desired acceleration.
            a_des = zeros(length(v_des), 1);
            a_des(1) = (v_des(2) - v_des(1)) / st;
            a_des(end) = (v_des(end) - v_des(end - 1)) / st;
            for i = 2:1:(length(a_des) - 1)
                a_des(i) = (v_des(i + 1) - v_des(i - 1)) / (2 * st);
            end

            % Calculate feedforward force depending on internal model
            if intModel == "full" || intModel == "no impedance"
                % Get desired ball angle trajectory (rad) and compute
                % feedforward force
                theta_des = int_mdl_output(:, 2);
                f = M * a_des - m * g * theta_des;
            elseif intModel == "rigid body"
                f = (m + M) * a_des;
            elseif intModel == "slow"
                f = M * a_des;
            elseif intModel == "fast"
                f = m * a_des;
            end
        else
            % If not using feedforward force, set f to zero vector
            x_des = xc';
            v_des = vc';
            f = zeros(length(vc), 1);
        end

        % Calculate control input to actual (external) system
        if intModel ~= "no impedance"
            u = f' + b * vc + k * xc;
        else
            u = f' + b * v_des' + k * x_des';
        end

        if ~impedance
            u = f';
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iv): Simulate system response
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~impedance
            ic = zeros(4,1);    % initial conditions
            [~,output] = ode45(@(t, y) nonlinSysEqnsNoImp(t,y,u,tc),tc,ic);

            % Assign outputs to variables
            pos_sim     = output(:, 1);
            theta_sim   = output(:, 2);
            vel_sim     = output(:, 3);
            omega_sim   = output(:, 4);
        else 
            % Reject profiles that end before the pendulum is released
            if pendIndex >= lentc
                startIndex = 0;
                output = 100000*ones(length(tc),6);
                vc = [];
                dur_corr = 0;
                invalidFlag = true;
                return
            end   
            
            % Divide time and input vectors into two sections, before and
            % after time that pendulum starts to move (when cart velocity
            % exceeds 0.1 m/s)
            t1 = tc(1:pendIndex);
            u1 = u(1:pendIndex);
            t2 = tc(pendIndex + 1:end) - tc(pendIndex + 1);
            u2 = u(pendIndex + 1:end);

            % Simulate rigid-body system for first portion of movement
            output1 = lsim(sysRigid, u1, t1);

            % Capture values at end of first section to use as initial
            % condition input to second section
            X0 = output1(end, :);

            % When pendulum is released, simulate full system for the
            % remainder of the move duration. Use values at end of previous
            % simulation as initial conditions for this simulation
            if simVersion == "linear"
                output2 = lsim(extSys, u2, t2, X0);
            elseif simVersion == "nonlinear"
                if impedance
                    [~,output2] = ode45(@(t,y) nonlinSysEqns(t,y,u2,t2,b,k),...
                        t2,X0);
                else
                    [~,output2] = ode45(@(t, y) nonlinSysEqnsNoImp(t,y,u2,t2),...
                        t2,X0);
                end
            end

            % Concatenate outputs and assign to variables
            output = [output1; output2];
            pos_sim     = output(:, 1);
            theta_sim   = output(:, 2);
            vel_sim     = output(:, 3);
            omega_sim   = output(:, 4);

        end

        % Take approximative derivative of cart velocity to get cart
        % acceleration profile. Take central difference at all interior
        % points, and forward and backward differences at first and last
        % points, respectively.
        acc_sim = zeros(length(vel_sim),1);
        acc_sim(1) = (vel_sim(2) - vel_sim(1))/st;
        acc_sim(end) = (vel_sim(end) - vel_sim(end-1))/st;
        for i = 2:1:length(acc_sim)-1
            acc_sim(i)  = (vel_sim(i+1) - vel_sim(i-1))/(2*st);
        end

        % Take approximate derivative of ball angular velocity in same way 
        % as cart to get ball angular acceleration profile.
        alpha_sim           = zeros(length(omega_sim),1);
        alpha_sim(1)        = (omega_sim(2) - omega_sim(1))/st;
        alpha_sim(end)      = (omega_sim(end) - omega_sim(end-1))/st;
        for i = 2:1:length(omega_sim)-1
            alpha_sim(i)    = (omega_sim(i+1) - omega_sim(i-1))/(2*st);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iv): Trim simulated profiles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Trim simulated output at first point where cart velocity exceeds 
        % experimental start velocity
        startIndex = find(vel_sim > vStart, 1);
        if ~isempty(startIndex)
            pos_sim     = pos_sim(startIndex:end);
            theta_sim   = theta_sim(startIndex:end);
            vel_sim     = vel_sim(startIndex:end);
            omega_sim   = omega_sim(startIndex:end);
            acc_sim     = acc_sim(startIndex:end);
            alpha_sim   = alpha_sim(startIndex:end);
            vc          = vc(startIndex:end);
        else
            startIndex  = 0;
        end

        % Find the last time the simulated velocity exceeds the threshold
        % and is decreasing. If it never exceeds the threshold, just take
        % the total length.
        corrVEnd = find((vel_sim > 0.1) & (acc_sim < 0), 1, 'last');
        if ~isempty(corrVEnd)
            vel_corr = vel_sim(1:corrVEnd);
        else
            vel_corr = vel_sim;
        end
        dur_corr = length(vel_corr)*st;

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
        end
        
        % Re-assign trimmed kinematic profiles to output array
        output = [pos_sim,theta_sim,vel_sim,omega_sim,acc_sim,alpha_sim];
        
    end

end
