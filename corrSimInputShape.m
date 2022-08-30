% Input Shaping Calculation and Simulation
%
% This function generates 2 or 4 input shaping impulses based on a system
% with 1 or 2 modes of oscillation, respectively. The impulses will be
% symmetric if the system is lossless (zeta = 0). If the system is
% underdamped (0 < zeta < 1), the impulses will be asymmetrical. If the
% system is critically or overdamped (zeta >= 1), input shaping will not work.

function [output, dur_corr] = corrSimInputShape(sys,Td,zeta,tDesSim,xEnd,vStart,st)
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (i): Generate Input Shaping impulses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate impulse times
    t1      = 0;                                                            
    t2      = Td/2;
    t2      = round(t2,mydigits(st));                                 % Round t2 to nearest time step
    tf      = tDesSim-t2;                                           % Duration of pre-convolved profile
    t       = 0:st:tf;                                                      
    tc      = 0:st:tDesSim;                                         % Time vector for convolved profile
    lentc   = length(tc);

    if tf <= 0
        output = 100000*ones(length(tc),6);
        dur_corr = NaN;
        return
    end

    % Generate impulse amplitudes. fa11 and fa12 are multiplicative
    % factors to reproduce human subject execution error.
    K       = exp(-zeta*pi/sqrt(1-zeta^2));     
    A1      = 1/(1+K);                                         % Amplitude of first impulse
    A2      = K/(1+K);                                         % Amplitude of second impulse  

    % Create vector for impulses
    pulses                  = zeros(1,length(t));                   % Initialize vector for impulses
    pulses(int16(1+t1/st))  = A1;                                   % Impulse 1 at t1
    pulses(int16(1+t2/st))  = A2;                                   % Impulse 2 at t2
    pulsetimes              = find(pulses);                         % Find indices of all pulses

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (ii): Shape desired input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Minimum-jerk velocity profile
    v = xEnd*(30*(t/tf).^2 - 60*(t/tf).^3 + 30*(t/tf).^4)/tf;

    % Convolve min jerk velocity profile with input shaping impulses
    vc = conv(v,pulses);    

    % Make time and velocity arrays same length
    if length(vc) > lentc
        vc = vc(1:lentc);                                               % Trim vector after desired time
    else
        vc(lentc) = 0;                                                  % Pad end of array with zeros
    end
    lenvc = length(vc);

    % Take approximative derivative of convolved velocity profile to
    % get input acceleration profile. Take central difference at all
    % interior points, forward difference at first point, and backward
    % difference at last point.
    ac         = zeros(1,lenvc);                                        % Initialize acceleration vector
    ac(1)      = (vc(2) - vc(1))/st;                                    % Forward difference
    ac(end)    = (vc(end) - vc(end-1))/st;                              % Backward difference

    for i = 2:1:length(ac)-1
        ac(i)  = (vc(i+1) - vc(i-1))/(2*st);                            % Central difference
    end

    % Make time and acceleration arrays same length
    if length(ac) > lentc
        ac = ac(1:lentc);                                               % Trim vector after desired time
    else
        ac(lentc) = 0;                                                  % Pad end of array with zeros
    end

    % Use desired acceleration trajectory as input
    u = ac;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (iii): Simulate system response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Simulate linearized system
    output  = lsim(sys,u,tc);

    % Assign outputs to variables
    pos_sim     = output(:,1);
    theta_sim   = output(:,2);
    vel_sim     = output(:,3);
    omega_sim   = output(:,4);

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

%     % Trim output at first point where cart velocity exceeds threshold
%     startIndex = find(vel_sim > vStart, 1);
%     if ~isempty(startIndex)
%         pos_sim     = pos_sim(startIndex:end);
%         theta_sim   = theta_sim(startIndex:end);
%         vel_sim     = vel_sim(startIndex:end);
%         omega_sim   = omega_sim(startIndex:end);
%         acc_sim     = acc_sim(startIndex:end);
%         alpha_sim   = alpha_sim(startIndex:end);
%     end
% 
%     % For duration-minimum velocity correlation, want to find what the
%     % duration of the simulated trial would be if it were trimmed the
%     % same way as experimental trials. Create "dummy" velocity variable
%     % and trim it for this purpose
%     corrVEnd = find((vel_sim > 0.1) & (acc_sim < 0), 1, 'last');
%     if ~isempty(corrVEnd)
%         vel_corr = vel_sim(1:corrVEnd);
%     else
%         vel_corr = vel_sim;
%     end
%     dur_corr = length(vel_corr)*st;

    dur_corr = length(vel_sim)*st;

    % Re-assign trimmed kinematic profiles to output array
    output = [pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim];

end
