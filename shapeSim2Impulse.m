% 2 input shaping impulses
%
% This function generates 2 input shaping impulses based on a system with
% one mode of oscillation. The impulses will be symmetric if the system is
% lossless (zeta = 0). If the system is underdamped (0 < zeta < 1), the
% impulses will be asymmetrical. If the system is critically or overdamped
% (zeta >= 1), input shaping will not work.

function output = shapeSim2Impulse(b,k,sys,Td1,Td2,zeta1,zeta2,tdes,xf,...
        vStart,fa11,fa12,fa21,fa22,p,q,r,s,st,shift,forwardF,linearized)
    
    % Generate impulse times
    t1      = 0;                                                            
    t2      = Td1/2;
    t2      = round(t2,digits(st));                                         % Round t2 to nearest time step
    tf      = tdes-t2;                                                      % Duration of pre-convolved profile
    t       = 0:st:tf;                                                      
    tc      = 0:st:tdes;                                                    % Time vector for convolved profile
    lentc   = length(tc);
    
    if tf <= 0
       output = 100000*ones(length(tc),6); 
       return
    end
    
    % Generate impulse amplitudes. a and b are multiplicative factors to
    % reproduce human subject execution error.
    K       = exp(-zeta1*pi/sqrt(1-zeta1^2));     
    A1      = fa11*1/(1+K);                                                  % Amplitude of first impulse
    A2      = fa12*K/(1+K);                                                  % Amplitude of second impulse  
    
    % Create vector for impulses
    pulses                  = zeros(1,length(t));                           % Initialize vector for impulses
    pulses(int16(1+t1/st))  = A1;                                           % Impulse 1 at t1
    pulses(int16(1+t2/st))  = A2;                                           % Impulse 2 at t2
    pulsetimes              = find(pulses);                                 % Find indices of all pulses
    
    % Impulse timing shifts:
    % Offset and shift p and q values so they each range from max to min
    % delay in terms of number of array elements. Wrap in int16 to make
    % sure indices are integer values.
    p = int16((p-2)*shift/st);
    q = int16((q-2)*shift/st); 
    
    % Keep original value if iteration is unshifted.
    % If p is positive, shift 1st pulse forward by shift time
    if p > 0
        % Identify index of shifted impulse.
        pulseindex = pulsetimes(1) + p;

        % Set amplitude at shifted pulse to value of original impulse
        pulses(pulseindex) = pulses(pulsetimes(1));

        % Set amplitude at original impulse time to zero
        pulses(pulsetimes(1)) = 0;
        
    % If p is negative, keep original value of 1st pulse and shift other 
    % pulse forward. 
    elseif p < 0
        
        % Impulse #2
        pulseindex              = pulsetimes(2) - p;
        pulses(pulseindex)      = pulses(pulsetimes(2));
        pulses(pulsetimes(2))   = 0;
        
        % Find indicies of all pulses again so that these become the new
        % indices for further shifts
        pulsetimes = find(pulses);
        
    end
    
    % Keep original value if iteration is unshifted.
    % Otherwise, shift forward or backward by shift time
    if q ~= 0
        % Identify index of shifted impulse. If shifted time is before
        % movementonset, set to t=0.
        pulseindex = pulsetimes(2) + q;
        if pulseindex <= 0
            pulseindex = 1;
        end
        
        % Set amplitude at shifted pulse to value of original impulse
        pulses(pulseindex) = pulses(pulsetimes(2));

        % Set amplitude at original impulse time to zero
        pulses(pulsetimes(2)) = 0;
    end

%%%%%%%%% Shape desired input and simulate system response %%%%%%%%%%%%%%%
    
    % Generate desired input using minimum-jerk profile (min jerk of hand)
    v = xf*(30*(t/tf).^2 - 60*(t/tf).^3 + 30*(t/tf).^4)/tf;

    % Convolve min. jerk velocity profile with input shaping impulses
    vc      = conv(v,pulses);
    lenvc   = length(vc);

    % Make time and velocity arrays same length
    if lenvc > lentc
        vc = vc(1:lentc);                                                   % Trim vector after desired time
    else
        vc(lentc) = 0;                                                      % Pad end of array with zeros
    end
    
    lenvc = length(vc);                                                     % Redefine length after trimming or padding

    if forwardF
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
        
    else
        % If no force input term, use desired vleocity trajectory as input
        u = vc;
    end

    % Simulate system response using convolved acceleration profile as
    % simulation input. Store results in output variable.
    if linearized
        % Simulate linearized system
        output  = lsim(sys,u,tc);
        
        % Assign outputs to variables
        pos_sim     = output(:,1);
        theta_sim   = output(:,2);
        vel_sim     = output(:,3);
        omega_sim   = output(:,4);
    else
        % Simulate system of nonlinear equations
        ic = zeros(6,1);
        [~,output] = ode45(@(t,y) nonlinSysEqns(t,y,u,tc,b,k), tc, ic);
        
        % Assign outputs to variables
        pos_sim     = output(:,2);
        theta_sim   = 360/(2*pi)*output(:,3);
        vel_sim     = output(:,5);
        omega_sim   = 360/(2*pi)*output(:,6);
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
    % to get ball acceleraiton profile.
    alpha_sim           = zeros(length(omega_sim),1);                       % Initialize ball acceleration vector
    alpha_sim(1)        = (omega_sim(2) - omega_sim(1))/st;                 % Forward difference
    alpha_sim(end)      = (omega_sim(end) - omega_sim(end-1))/st;           % Backward difference
    for i = 2:1:length(omega_sim)-1
        alpha_sim(i)    = (omega_sim(i+1) - omega_sim(i-1))/(2*st);         % Central difference
    end
    
    % Trim output at point where cart velocity reaches same start velocity
    % as experimental data
    simVStart     = find(vel_sim > vStart);
    if ~isempty(simVStart)
        startIndex  = simVStart(1);
        pos_sim     = pos_sim(startIndex:end);
        theta_sim   = theta_sim(startIndex:end);
        vel_sim     = vel_sim(startIndex:end);
        omega_sim   = omega_sim(startIndex:end);
        acc_sim     = acc_sim(startIndex:end);
        alpha_sim   = alpha_sim(startIndex:end);
    end
    
    % Re-assign trimmed kinematic profiles to "output" array
    output = [pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim];

end
