% Trims cup task experimental data from HapticMaster given particular cart
% position and velocity characteristics.
%
% Last updated: 2020-10-15    
    
clear all;

fileStr = 'experimental data/Data Collection 7-17-2019/CM/Block 4/CM_17Jul2019_16-47-55_trial_30.mat';
load(fileStr,'pos','theta','vel','omega','acc','alpha','t')

st          = .001;
delayMax    = 0.100;
start       = 276;
stop        = 371;

% Convert ball units to degrees
theta = theta*360/(2*pi);
omega = omega*360/(2*pi);
alpha = alpha*360/(2*pi);

% Offset cart position so it starts at x = 0
pos     = pos + 0.125;

% Trim vectors to specified start and stop indices
t       = t(start:stop);
pos     = pos(start:stop);
theta   = theta(start:stop);
vel     = vel(start:stop);
omega   = omega(start:stop);
acc     = acc(start:stop);
alpha   = alpha(start:stop);

% Shift experimental time vector and create linear time vector
t       = t - t(1);                                                     % Start time vector at t=0
tdes    = t(end) - t(1);                                                % Total experimental movement time duration
tc      = 0:st:tdes;

% Interpolate data to match simulated sample points
pos     = interp1(t,pos,tc);
theta   = interp1(t,theta,tc);
vel     = interp1(t,vel,tc);
omega   = interp1(t,omega,tc);
acc     = interp1(t,acc,tc);
alpha   = interp1(t,alpha,tc);

% Output interpolated time vector and kinematic data as column vectors
t = tc';
pos = pos';
theta = theta';
vel = vel';
omega = omega';
acc = acc';
alpha = alpha';

save('CM_4_31_Experimental.mat','t','pos','theta','vel','omega','acc','alpha');

figure();
plot(t,vel)