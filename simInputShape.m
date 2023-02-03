% Input Shaping Calculation and Simulation
%
% This function generates 2 or 4 input shaping impulses based on a system
% with 1 or 2 modes of oscillation, respectively. The impulses will be
% symmetric if the system is lossless (zeta = 0). If the system is
% underdamped (0 < zeta < 1), the impulses will be asymmetrical. If the
% system is critically or overdamped (zeta >= 1), input shaping will not work.

function [output, vcOut, sub1End, sub2Start, dur_corr] = simInputShape(b,k,ver,sys,...
        sysRigid,Td1,Td2,zeta1,zeta2,tExp,tDesSim,xEnd,vStart,vEnd,fa11,...
        fa12,fa21,fa22,p,q,r,s,st,shift,forwardF,simVersion,modes,...
        pendIndex,fitMethod)

    G = globalData();
    m = G.m;
    M = G.M;
    l = G.l;
    g = G.g;

    % DEBUG
%     disp("Inside simInputShape")
    
    % Reset flag for invalid trials
    flag = false;
    
    % Run first-pass simulation
    [output, startIndex, vcOut, sub1End, sub2Start, dur_corr] = sim();
    
    if flag
        % If first run produces a rejection, take the output and return to
        % the optimization script without running a second time
        return
    end
             
    % Add amount of time that was trimmed to pendulum release time
    pendIndex = pendIndex + startIndex;
    
    % Simulate again with lengthened duration
    [output, ~, vcOut, sub1End, sub2Start, dur_corr] = sim();
    
    function [output, startIndex, vcOut, sub1End, sub2Start, dur_corr] = sim()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (i): Generate Input Shaping impulses
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % 1 mode (2 input shaping impulses)
        if modes == 1
            % Generate impulse times
            t1      = 0;                                                            
            t2      = Td1/2;
            t2      = round(t2,mydigits(st));                                 % Round t2 to nearest time step
            tf      = tDesSim-t2;                                           % Duration of pre-convolved profile
            t       = 0:st:tf;                                              % Time vector: pre-convolution
            tc      = 0:st:tDesSim;                                         % Time vector: post-convolution
            lentc   = length(tc);

            if tf <= 0
                startIndex = 0;
                output = 100000*ones(length(tc),6);
                vc = [];
                vcOut = [];
                sub1End = 99;
                sub2Start = 99;
                dur_corr = 0;
                flag = true;
                return
            end

            % Generate impulse amplitudes. fa11 and fa12 are multiplicative
            % factors to reproduce human subject execution error.
            K       = exp(-zeta1*pi/sqrt(1-zeta1^2));     
            A1      = fa11*1/(1+K);                                         % Amplitude of first impulse
            A2      = fa12*K/(1+K);                                         % Amplitude of second impulse  
            
            % Create vector for impulses
            pulses                  = zeros(1,length(t));                   % Initialize vector for impulses
            pulses(int16(1+t1/st))  = A1;                                   % Impulse 1 at t1
            pulses(int16(1+t2/st))  = A2;                                   % Impulse 2 at t2
            pulsetimes              = find(pulses);                         % Find indices of all pulses

        % 2 modes (4 input shaping impulses)
        else
            % Generate impulse times
            t1      = 0;                                                        % Set t1 to t=0 for both modes
            t21     = Td1/2;                                                    % Second impulse for first mode
            t22     = Td2/2;                                                    % Second impulse for second mode
            toffset = t21+t22;
            tf      = tDesSim-toffset;                                             % Duration of pre-convolved profile
            t       = 0:st:tf;
            tc      = 0:st:tDesSim;                                                % Time vector for convolved profile
            lentc   = length(tc);

%             % DEBUG
%             tDesSim
%             toffset
%             t21
%             t22
%             lent = length(t)
%             lentc
%             Td1
%             Td2
%             zeta1
%             zeta2
%             tDesSim
%             xEnd

            if tf <= 0
                startIndex = 0;
                output = 100000*ones(length(tc),6);
                vc = [];
                vcOut = [];
                sub1End = 99;
                sub2Start = 99;
                dur_corr = 0;
                flag = true;
                return
            end

            % Generate first mode impulse amplitudes
            K1      = exp(-zeta1*pi/sqrt(1-zeta1^2));     
            A11     = fa11*1/(1+K1);                                            % Amplitude of first impulse
            A21     = fa12*K1/(1+K1);                                           % Amplitude of second impulse

            % Generate second mode impulse amplitudes
            K2      = exp(-zeta2*pi/sqrt(1-zeta2^2));     
            A12     = fa21*1/(1+K2);                        
            A22     = fa22*K2/(1+K2);             

            % Create vector for first mode impulses
            t21     = round(t21,mydigits(st));                                    % Round t21 to nearest time step
            imp1    = zeros(1,length(t));                                       % Initialize imp vector for input impulses
            imp1(int16(1+t1/st))    = A11;                                      % Impulse 1 at t1
            imp1(int16(1+t21/st))   = A21;                                      % Impulse 2 at t21

            % Create vector for second mode impulses
            t22     = round(t22,mydigits(st));                                    % Round t22 to nearest time step
            imp2    = zeros(1,length(t));                                       % Initialize imp vector for input impulses
            imp2(int16(1+t1/st))    = A12;                                      % Impulse 1 at t1
            imp2(int16(1+t22/st))   = A22;                                      % Impulse 2 at t22

            % Convolve two sets of impulses
            pulses  = conv(imp1,imp2);

            % Find indicies of all pulses
            pulsetimes = find(pulses);

            % Trim trailing zeros from vector of pulses
            pulses = pulses(1:pulsetimes(end));

%             % DEBUG
%             pulsetimes
%             length(pulses)

            if forwardF
                % Create dynamic system of internal model
                [internal_sys,~,~,~,~,~,~,~,~] = sysCreate(b,k,forwardF,ver,true,false);

                % DEBUG
%                 disp("Back to simInputShape")
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (ii): Shape desired input
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Minimum-jerk profiles
%         x = xEnd * (10 * (t / tf).^3 - 15 * (t / tf).^4 + 6 * (t / tf).^5);
        v = xEnd * (30 * (t / tf).^2 - 60 * (t / tf).^3 + 30 * (t / tf).^4) / tf;
        a = xEnd * (60 * (t / tf) - 180 * (t / tf).^2 + 120 * (t / tf).^3) / (tf ^ 2);

        % Convolve min jerk profiles with input shaping impulses
%         xc = conv(x, pulses);
        vc = conv(v, pulses);
        ac = conv(a, pulses);

        % DELETE this later if not needed %%
        % Populate array to hold individual velocity submovement data
        vcOut = zeros(length(pulsetimes), length(vc));
        for pulseCount = 1:length(pulsetimes)
            tempPulses = zeros(1,length(pulses));
            currentPulse = pulsetimes(pulseCount);
            tempPulses(currentPulse) = pulses(currentPulse);
            vcOut(pulseCount,:) = conv(v,tempPulses);
        end
        
        % Identify timing of submovement 1 end and submovement 2 start
        if length(pulsetimes) == 2
            sub1End = st*find(vcOut(1,:), 1, 'last');
            sub2Start = st*find(vcOut(2,:), 1, 'first');
        elseif length(pulsetimes) == 4
            sub1End = st*find(vcOut(1,:), 1, 'last');
            sub2Start = st*find(vcOut(3,:), 1, 'first');
        else
            lenPulses = length(pulsetimes);
            sub1End = 99;
            sub2Start = 99;
            disp('sub1End & sub2Start not assigned!')
            disp(' ')
        end
        % END DELETE %%

        % Make time and kinematics arrays same length, if rounding error
        if length(vc) > lentc
            % Trim vectors after desired time
            vc = vc(1:lentc);
            ac = ac(1:lentc);
        else
            % Pad end of array with zeros
            vc(lentc) = 0;
            ac(lentc) = 0;
        end
        lenvc = length(vc);

%         % DEBUG
%         disp("Length of tc: ")
%         length(tc)
%         disp("Length of vc: ")
%         length(vc)

        xc = cumtrapz(tc, vc);
        xc(end) = xc(end - 1); % correct cumtrapz calculating 0 for last index

%         % Take approximative derivative of convolved velocity profile to
%         % get input acceleration profile. Take central difference at all
%         % interior points, forward difference at first point, and backward
%         % difference at last point.
%         ac         = zeros(1,lenvc);                                        % Initialize acceleration vector
%         ac(1)      = (vc(2) - vc(1))/st;                                    % Forward difference
%         ac(end)    = (vc(end) - vc(end-1))/st;                              % Backward difference
% 
%         for i = 2:1:length(ac)-1
%             ac(i)  = (vc(i+1) - vc(i-1))/(2*st);                            % Central difference
%         end
% 
%         % Make time and acceleration arrays same length
%         if length(ac) > lentc
%             ac = ac(1:lentc);                                               % Trim vector after desired time
%         else
%             ac(lentc) = 0;                                                  % Pad end of array with zeros
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iii): Simulate feedforward system response
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if modes == 2
            % DEBUG
%             disp('Simulating feedforward response!')
%             disp("Size of xc: ")
%             disp(size(xc))
%             disp("Size of vc: ")
%             disp(size(vc))

            % Simulate response of internal model when given shaped input
            u = (b / k * vc) + xc ;
            int_mdl_output = lsim(internal_sys, u, tc);

            % DEBUG
%             figure();
%             plot(tc, u)
%             title("Internal model input")
% 
%             figure();
%             plot(tc, xc)
%             title("Convolved position")
% 
%             figure();
%             plot(tc, vc)
%             title("Convolved velocity")
% 
%             figure();
%             plot(tc, int_mdl_output(:,1))
%             title("Internal simulated cart position")
% 
%             figure();
%             plot(tc, int_mdl_output(:,2))
%             title("Internal simulated ball angle")
% 
%             figure();
%             plot(tc, int_mdl_output(:,3))
%             title("Internal simulated cart velocity")
% 
%             figure();
%             plot(tc, int_mdl_output(:,4))
%             title("Internal simulated ball angular velocity")
% 
%             disp("Internal system: ")
%             internal_sys

            theta_des = deg2rad(int_mdl_output(:, 2));
            v_des = int_mdl_output(:, 3);

            % DEBUG
%             disp("size of theta_des:")
%             size(theta_des)
%             disp("size of v_des:")
%             size(v_des)

            % Differentiate desired velocity to get desired acceleration 
            a_des = zeros(length(v_des), 1);                                % Initialize acceleration row vector
            a_des(1) = (v_des(2) - v_des(1)) / st;                          % Forward difference
            a_des(end) = (v_des(end) - v_des(end - 1)) / st;                % Backward difference
            for i = 2:1:(length(a_des) - 1)
                a_des(i) = (v_des(i + 1) - v_des(i - 1)) / (2 * st);            % Central difference
            end
    
            % Calculate feedforward force
            f = M * a_des - m * g * theta_des;
%             f = M * a_des;
%             f = M * ac' - m * g * theta_des;
%             f = M * ac';
    
            % Calculate system control input including feedforward force
            u = (1 / k * f') + (b / k * vc) + xc ;

            % DEBUG
%             disp("Input w/ feedforward force dims: ")
%             size(u)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % (iv): Simulate system response
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if simVersion == "linear"
            % Simulate linearized system
            output  = lsim(sys, u, tc);
 
            % Assign outputs to variables
            pos_sim     = output(:, 1);
            theta_sim   = output(:, 2);
            vel_sim     = output(:, 3);
            omega_sim   = output(:, 4);

        % Simulate system with pendlulum lock bug
        elseif simVersion == "pendLock"
            
            % Reject profiles that end before the pendulum is released
            if pendIndex >= lentc
                startIndex = 0;
                output = 100000*ones(length(tc),6);
                vc = [];
                vcOut = [];
                sub1End = 99;
                sub2Start = 99;
                dur_corr = 0;
                flag = true;
                return
            end   
            
            % Divide time and input vectors into two sections, before and
            % after time that pendulum starts to move (when cart velocity
            % exceeds 0.1 m/s)
            t1 = tc(1:pendIndex);
            u1 = u(1:pendIndex);
            t2 = tc(pendIndex+1:end)-tc(pendIndex+1);                       % Shift value so first t = 0
            u2 = u(pendIndex+1:end);

            % DEBUG
%             disp("sizes of t1, u1, t2, u2: ")
%             size(t1)
%             size(u1)
%             size(t2)
%             size(u2)

            % Simulate rigid-body system for first portion of movement
            output1 = lsim(sysRigid, u1, t1);

            % Capture values at end of first section to use as initial
            % condition input to second section
            if forwardF
%                 if length(output1(1, :)) ~= 4
% 
%                     X0 = [output1(end,5), output1(end,1),...                % xdes, x,
%                         output1(end,2)*2*pi/360, vc(pendIndex), ...         % theta (rad), xdes_dot = vc,
%                         output1(end,3), output1(end,4)*2*pi/360];           % x_dot, theta_dot (rad/s)
%                 else
                    X0 = output1(end, :);
%                     X0 = [output1(end, 1);                                  % x
%                           output1(end, 2);                                  % theta (rad)
%                           output1(end, 3);                                  % x_dot
%                           output1(end, 4)];                                 % theta_dot (rad/s)
%                 end
            else
%                 X0 = [output1(end,5), output1(end,1),...                    % xdes, x,
%                     output1(end,2)*2*pi/360, vc(pendIndex),...              % theta (rad), xdes_dot = vc,
%                     output1(end,3), output1(end,4)*2*pi/360];               % x_dot, theta_dot (rad/s)
            end

            % When pendulum is released, simulate full system for the
            % remainder of the move duration. Use values at end of previous
            % simulation as initial conditions for this simulation
            output2 = lsim(sys, u2, t2, X0);

            % Concatenate outputs
            output = [output1;
                      output2];

            % Assign outputs to variables
            pos_sim     = output(:, 1);
            theta_sim   = output(:, 2);
            vel_sim     = output(:, 3);
            omega_sim   = output(:, 4);

        elseif simVersion == "nonlinear"
            % Simulate system of nonlinear equations
            ic = zeros(6,1);
            [~,output] = ode45(@(t,y) nonlinSysEqns(t,y,u,tc,b,k), tc, ic);

            % Assign outputs to variables
            pos_sim     = output(:,2);
            theta_sim   = output(:,3);
            vel_sim     = output(:,5);
            omega_sim   = output(:,6);
        end

        % Take approximative derivative of cart velocity to get cart
        % acceleration profile. Take central difference at all interior points,
        % forward difference at first point, and backward difference at last
        % point.
        acc_sim         = zeros(length(vel_sim),1);                         % Initialize cart acceleration vector
        acc_sim(1)      = (vel_sim(2) - vel_sim(1))/st;                     % Forward difference
        acc_sim(end)    = (vel_sim(end) - vel_sim(end-1))/st;               % Backward difference
        for i = 2:1:length(acc_sim)-1
            acc_sim(i)  = (vel_sim(i+1) - vel_sim(i-1))/(2*st);             % Central difference
        end

        % Take approximate derivative of ball angular velocity in same way 
        % as cart to get ball angular acceleration profile.
        alpha_sim           = zeros(length(omega_sim),1);                   % Initialize ball acceleration vector
        alpha_sim(1)        = (omega_sim(2) - omega_sim(1))/st;             % Forward difference
        alpha_sim(end)      = (omega_sim(end) - omega_sim(end-1))/st;       % Backward difference
        for i = 2:1:length(omega_sim)-1
            alpha_sim(i)    = (omega_sim(i+1) - omega_sim(i-1))/(2*st);     % Central difference
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
            vcOut       = vcOut(:, startIndex:end);
        else
            startIndex  = 0;
        end

        % For duration-minimum velocity correlation, want to find what the
        % duration of the simulated trial would be if it were trimmed the
        % same way as experimental trials. Create "dummy" velocity variable
        % and trim it for this purpose
        lenExp = length(0:st:tExp);
        lenSim = length(vel_sim);
%         % If simulated profile is shorter than experimental profile, take
%         % duration as full simulated length
%         if lenSim <= lenExp
%             dur_corr = lenSim*st;
%         % If simulated velocity at the end of the experimental trial
%         % duration is already below the threshold, trim it there
%         elseif vel_sim(lenExp) < 0.1
%             vel_corr = vel_sim(1:lenExp);
%             dur_corr = length(vel_corr)*st;
        % Otherwise, find the last time the simulated velocity exceeds
        % the threshold and is decreasing. If it never exceeds the
        % threshold, just take the
        % total length.
%         else
            corrVEnd = find((vel_sim > 0.1) & (acc_sim < 0), 1, 'last');
            if ~isempty(corrVEnd)
                vel_corr = vel_sim(1:corrVEnd);
            else
                vel_corr = vel_sim;
            end
            dur_corr = length(vel_corr)*st;
%         end
        
        % Trim simulated output at last point where cart velocity exceeds 
        % experimental data end velocity 
%         simVEnd = find(vel_sim > vEnd, 1, 'last');
%         if ~isempty(simVEnd)
%             pos_sim     = pos_sim(1:simVEnd);
%             theta_sim   = theta_sim(1:simVEnd);
%             vel_sim     = vel_sim(1:simVEnd);
%             omega_sim   = omega_sim(1:simVEnd);
%             acc_sim     = acc_sim(1:simVEnd);
%             alpha_sim   = alpha_sim(1:simVEnd);
%         end

        if fitMethod == "fitToAverage"
            % When fitting to the average of trials, trim end of trial at
            % velocity-based index, not to match length of individual trial
            simVEnd = find(vel_sim > 0.02, 1, 'last');
            if ~isempty(simVEnd)
                pos_sim     = pos_sim(1:simVEnd);
                theta_sim   = theta_sim(1:simVEnd);
                vel_sim     = vel_sim(1:simVEnd);
                omega_sim   = omega_sim(1:simVEnd);
                acc_sim     = acc_sim(1:simVEnd);
                alpha_sim   = alpha_sim(1:simVEnd);
            end
            
        else
            % If simulated profile is longer than experimental profile,
            % trim to same length
            lenExp = length(0:st:tExp);
            if length(vel_sim) > lenExp
                pos_sim     = pos_sim(1:lenExp);
                theta_sim   = theta_sim(1:lenExp);
                vel_sim     = vel_sim(1:lenExp);
                omega_sim   = omega_sim(1:lenExp);
                acc_sim     = acc_sim(1:lenExp);
                alpha_sim   = alpha_sim(1:lenExp);
            end
            if max(size(vcOut)) > lenExp
                vcOut       = vcOut(:,1:lenExp);
            end
        end
        
        % Re-assign trimmed kinematic profiles to output array
        output = [pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim];
        
    end

end
