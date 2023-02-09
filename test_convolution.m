% Convolve pulses with minimum jerk profiles. View any differences from
% non-commutativity of convolution.

Td1 = 1.4850;
Td2 = 0.5227;
zeta1 = 0.0061;
zeta2 = 0.1600;
tDesSim = 1.5324;
st = 0.001;

t1      = 0;                                                        % Set t1 to t=0 for both modes
t21     = Td1/2;                                                    % Second impulse for first mode
t22     = Td2/2;                                                    % Second impulse for second mode
toffset = t21+t22;
tf      = tDesSim-toffset;                                             % Duration of pre-convolved profile
t       = 0:st:tf;
tc      = 0:st:tDesSim;                                                % Time vector for convolved profile
lentc   = length(tc);

% Generate first mode impulse amplitudes
K1      = exp(-zeta1*pi/sqrt(1-zeta1^2));     
A11     = 1/(1+K1);                                            % Amplitude of first impulse
A21     = K1/(1+K1);                                           % Amplitude of second impulse

% Generate second mode impulse amplitudes
K2      = exp(-zeta2*pi/sqrt(1-zeta2^2));     
A12     = 1/(1+K2);                        
A22     = K2/(1+K2);             

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

xEnd = 0.24;

% Minimum-jerk profiles
x = xEnd * (10 * (t / tf).^3 - 15 * (t / tf).^4 + 6 * (t / tf).^5);
v = xEnd * (30 * (t / tf).^2 - 60 * (t / tf).^3 + 30 * (t / tf).^4) / tf;
a = xEnd * (60 * (t / tf) - 180 * (t / tf).^2 + 120 * (t / tf).^3) / (tf ^ 2);

% Direct convolvution
xc = conv(x, pulses);
vc = conv(v, pulses);
ac = conv(a, pulses);

% Differentation after convolution
% Take approximative derivative of convolved velocity profile to
% get input acceleration profile. Take central difference at all
% interior points, forward difference at first point, and backward
% difference at last point.
ad         = zeros(1,length(vc));                                        % Initialize adceleration vector
ad(1)      = (vc(2) - vc(1))/st;                                    % Forward difference
ad(end)    = (vc(end) - vc(end-1))/st;                              % Backward difference

for i = 2:1:length(ad)-1
    ad(i)  = (vc(i+1) - vc(i-1))/(2*st);                            % Central difference
end

% Take integral of vc to get cumulative position
xd = cumtrapz(tc, vc);

% Plot both accelerations to compare
% figure();
% plot(tc,ac)
% hold on;
% plot(tc,ad)

% figure();
% title("Minimum jerk acceleration")
% plot(t,a)

figure();
plot(t,x)
title("Minimum jerk position")

figure();
plot(t,x)
title("Minimum jerk position")

figure();
plot(tc, xc);
title("Convolved position xc")

figure();
plot(tc, xd);
title("Integrated position xd")

