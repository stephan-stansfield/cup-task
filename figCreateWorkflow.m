% Produce figures to help explain input shaping simulation and fitting
% process for work flow diagram.

close all;
clear variables;
clc;

%% Show "raw" experimental profile

close all;

% Define global variables
global m M l g;
m           = 1.1;                                                          % Mass of ball (kg)
M           = 1.9;                                                          % Mass of cup (kg)
l           = 0.5;                                                          % Length of the massless pendulum rod (m)
g           = 9.8;                                                          % Gravity (m/s^2)

% Choose experimental trial
% num         = 42;
num         = 45;
subjStr     = "NO";
trialDate   = "_15Jul2019";
trialStr    = "_13-17-36_trial_";
blockStr    = "NOb3.mat";
numStr      = num2str(num);

% Construct file name and load file of one experimental trial
fileStr = strcat(subjStr,trialDate,trialStr,numStr);
load(fileStr,'pos','theta','vel','omega','acc','alpha','t')

% Convert ball units from radians to degrees
theta = rad2deg(theta);
omega = rad2deg(omega);
alpha = rad2deg(alpha);

% Create figures of raw kinematic variables: cart position and velocity
figure();
plot(t,pos,'LineWidth',2)
set(gca,'FontSize',16)
title('Untrimmed Experimental Cart Position Trajectory','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Cart Position (m)','FontSize',18)

figure();
plot(t,vel,'LineWidth',2)
set(gca,'FontSize',16)
title('Untrimmed Experimental Cart Velocity Profile','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Cart Velocity (m/s)','FontSize',18)

%%

% % Plot other untrimmed kinematic variables
% figure();
% plot(t,acc,'LineWidth',2)
% set(gca,'FontSize',16)
% title('Untrimmed Experimental Cart Acceleration Profile','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Cart Acceleration (m/s^2)','FontSize',18)
% 
% figure();
% plot(t,theta,'LineWidth',2)
% set(gca,'FontSize',16)
% title('Untrimmed Experimental Ball Angular Trajectory','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Ball Angle (deg)','FontSize',18)
% 
% figure();
% plot(t,omega,'LineWidth',2)
% set(gca,'FontSize',16)
% title('Untrimmed Experimental Ball Angular Velocity Profile','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Ball Angular Velocity (deg/s)','FontSize',18)
% 
% figure();
% plot(t,alpha,'LineWidth',2)
% set(gca,'FontSize',16)
% title('Untrimmed Experimental Ball Angular Acceleration Profile','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Ball Angular Acceleration (deg/s^2)','FontSize',18)

%% Show trimmed experimental profile

% Load start and stop indices of corresponding trial. Note that trials
% are indexed from 0, so add 1 to access correct row in array.
load(blockStr,'Expression1')
start   = Expression1(num+1,1);
stop    = Expression1(num+1,2);

% Set time step
st = 0.001;

% Trim data
theta = deg2rad(theta);
omega = deg2rad(omega);
alpha = deg2rad(alpha);
[pos,theta,vel,omega,acc,alpha,tdes] = trimData(pos,theta,vel,omega,acc,alpha,t,st,start,stop);
xEnd = pos(end);
vEnd = vel(end);

t = 0:st:st*(length(pos)-1);

%{
% Plot trimmed kinematic variable profiles: cart position and velocity
figure();
plot(t,pos,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Cart Position Trajectory','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Cart Position (m)','FontSize',18)

figure();
plot(t,vel,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Cart Velocity Profile','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Cart Velocity (m/s)','FontSize',18)

%%

% Plot other trimmed kinematic variables
figure();
plot(t,acc,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Cart Acceleration Profile','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Cart Acceleration (m/s^2)','FontSize',18)

figure();
plot(t,theta,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Ball Angular Trajectory','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Ball Angle (deg)','FontSize',18)

figure();
plot(t,omega,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Ball Angular Velocity Profile','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Ball Angular Velocity (deg/s)','FontSize',18)

figure();
plot(t,alpha,'LineWidth',2)
set(gca,'FontSize',16)
title('Trimmed Experimental Ball Angular Acceleration Profile','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Ball Angular Acceleration (deg/s^2)','FontSize',18)
%}

%% Generate system and input shaping impulses - 2 impulses

% Simplified system with 2 impulses using rigid body internal model

% For trial 43 slow mode best-fit:
% b = 8.23;
% k = 250;
% ver = "slow";

% For trial 38 slow mode best-fit:
% b = 27.78;
% k = 518.5185195;
% ver = "slow";
% b = 30;
% k = 150;
% ver = "slow";

% For rigid body best-fit:
b = 7.71;
k = 106.524;
ver = "rigid body";

forwardF = true;
impedance = true;
print = false;
% [sys,Td,zeta,overdamped] = sysCreateSimple(b,k,ver,forwardF);

[sys,sysRigid,Td,Td1,Td2,zeta,zeta1,zeta2,overdamped] = sysCreate(b,k,...
        forwardF,ver,impedance,print);
    
% DEBUG
disp('Just got out of sysCreate')
Td
% Td1
% Td2
zeta
% zeta1
% zeta2

% Generate impulse times
t1      = 0;                                                            
t2      = Td/2;
t2      = round(t2,digits(st));                                         
tf      = tdes-t2;                                                      
t       = 0:st:tf;                                                      
tc      = 0:st:tdes;                                                    
lentc   = length(tc);

% Generate impulse amplitudes. a and b are multiplicative factors to
% reproduce human subject execution error.
fa1     = 1;
fa2     = 1;
K       = exp(-zeta*pi/sqrt(1-zeta^2));     
A1      = fa1*1/(1+K);                                                  
A2      = fa2*K/(1+K);                                                  

% Create vector for impulses
pulses                  = zeros(1,length(t));                           
pulses(int16(1+t1/st))  = A1;                                           
pulses(int16(1+t2/st))  = A2;                                           
pulsetimes              = st*find(pulses);
pulsePlot               = [A1, A2];

% Plot input shaping impulses
figure();
stem(pulsetimes,pulsePlot,'Marker','none','LineWidth',4)
set(gca,'FontSize',14)
title('Input Shaping Impulses - Rigid Body Internal Model','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Impulse Amplitude','FontSize',18)
xlim([0 1])
ylim([0 1])


%% Generate minimum jerk submovement shape

v = xEnd*(30*(t/tf).^2 - 60*(t/tf).^3 + 30*(t/tf).^4)/tf;

% Plot minimum jerk submovement velocity profile
figure();
plot(t,v,'LineWidth',10,'Color','#77AC30')
% set(gca,'FontSize',16)
% title('Minimum Jerk Velocity Submovement','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Velocity (m/s)','FontSize',18)
% set(gca,'xtick',[], 'ytick', [])
set(gca,'Visible','off')
box off

%% Convolve impulses with minimum jerk submovement

% Show minimum jerk velocity profile convolved with each individual impulse
% First impulse
pulse1 = pulses;
pulse1(pulsetimes(2)/st) = 0;
v1 = conv(v,pulse1);
if length(v1) > lentc
        v1 = v1(1:lentc);                                                   % Trim vector after desired time
end

figure();
plot(tc,v1,'LineWidth',2,'Color','#0072BD')                                 % Blue
set(gca,'FontSize',16)
hold on;
stem(pulsetimes(1),pulsePlot(1),'Marker','none','LineWidth',4,'Color','#0072BD')
plot(t,v,'LineWidth',2,'Color','#77AC30')                                   % Green: Original profile
hold off;
title('Convolution of First Impulse with Submovement','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylim([0 1])

% Second impulse
pulse2 = pulses;
pulse2(pulsetimes(1)/st) = 0;
v2 = conv(v,pulse2);
if length(v2) > lentc
        v2 = v2(1:lentc);                                                   % Trim vector after desired time
end

figure();
plot(tc,v2,'LineWidth',2,'Color','#A2142F')                                 % Red
set(gca,'FontSize',16)
hold on;
stem(pulsetimes(2),pulsePlot(2),'Marker','none','LineWidth',4,'Color','#A2142F')
plot(t,v,'LineWidth',2,'Color','#77AC30')                                   % Green: Original profile
hold off;
title('Convolution of Second Impulse with Submovement','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylim([0 1])

% Convolve minimum jerk profile with both impulses
vc      = conv(v,pulses);
if length(vc) > lentc
        vc = vc(1:lentc);                                                   % Trim vector after desired time
end
lenvc   = length(vc);

% Plot result of convolution between impulses and submovement
% (zero force velocity profile)
figure();
plot(tc,vc,'LineWidth',10,'Color','#7E2F8E')                                 % Purple                                
% set(gca,'FontSize',16)
% title('Zero-Force Hand Velocity Profile','FontSize',20)
% xlabel('Time (s)','FontSize',18)
% ylabel('Velocity (m/s)','FontSize',18)
set(gca,'Visible','off')
box off


% Plot result showing constituent impulses
figure();
plot(tc,vc,'LineWidth',2,'Color','#7E2F8E')                                 % Purple
hold on;
% stem(pulsetimes,pulsePlot,'Marker','none','LineWidth',4)
plot(tc,v1,'LineWidth',2,'Color','#0072BD')
plot(tc,v2,'LineWidth',2,'Color','#A2142F')
% stem(pulsetimes(1),pulsePlot(1),'Marker','none','LineWidth',4,'Color','#0072BD')
% stem(pulsetimes(2),pulsePlot(2),'Marker','none','LineWidth',4,'Color','#A2142F')
hold off;
set(gca,'FontSize',16)
title({'Summed Submovements:', 'Zero-Force Hand Velocity Profile'},'FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Velocity (m/s)','FontSize',18)
% ylim([0 1])

%% Simulate system using convolved profile as input

ac         = zeros(1,lenvc);                                            % Initialize acceleration vector
ac(1)      = (vc(2) - vc(1))/st;                                        % Forward difference
ac(end)    = (vc(end) - vc(end-1))/st;                                  % Backward difference

for i = 2:1:length(ac)-1
    ac(i)  = (vc(i+1) - vc(i-1))/(2*st);                                % Central difference
end

% Make time and acceleration arrays same length
if length(ac) > lentc
    ac = ac(1:lentc);                                                   % Trim vector after desired time
else
    ac(lentc) = 0;                                                      % Pad end of array with zeros
end

% Use desired acceleration trajectory as input
u = ac;

% Simulate system response using convolved acceleration profile as
% simulation input. Store results in output variable.
output      = lsim(sys,u,tc);

% Assign outputs to variables
pos_sim     = output(:,1);
theta_sim   = output(:,2);
vel_sim     = output(:,3);
omega_sim   = output(:,4);

% Plot output cart velocity
figure();
plot(tc,vel_sim,'LineWidth',2,'Color','#D95319')                            % Orange
set(gca,'FontSize',16)
title('Output Cart Velocity','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Velocity (m/s)','FontSize',18)

% Plot output cart position
figure();
plot(tc,pos_sim,'LineWidth',2)
set(gca,'FontSize',16)
title('Output Cart Position','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Position (m)','FontSize',18)

% Compute trimmed velocity profile
vStart = vel(1);
minA1 = 1;
minA2 = 1;
minP = 0;
minQ = 0;
tDesMax = 0;
shift = 0;

% For trial 43, slow mode w/ timing and amplitude variability best-fit:
% minP = -0.00164609053497942;
% minQ = -0.120164609;
% tdes = tdes + 0.099999993;
% minA1 = 0.988721022;
% minA2 = 0.921124724;
% tDesMax = 0;
% shift = 0.200;

% For trial 38, slow mode w/ timing and amplitude variability best-fit:
minP = 0.001646091;
minQ = -0.014814815;
tdessim = tdes + 0.099996611;
minA1 = 1.06975493;
minA2 = 0.903172648;
tDesMax = 0;
shift = 0.200;

% output = shapeSim2Impulse(sys, Td, zeta, tdes, xf, vStart, minA1,...
%                     minA2,minP, minQ, st, tDesMax, shift, forwardF);
%                 
% output = shapeSim2Impulse(b,k,sys,Td1,Td2,zeta1,zeta2,tdes,xf,...
%         vStart,fa11,fa12,fa21,fa22,p,q,r,s,st,shift,forwardF,linearized)

% Try with no variability
fa11 = 1;
fa12 = 1;
fa21 = 1;
fa22 = 1;
p = 0;
q = 0;
r = 0;
s = 0;
% p = 2;
% q = 2;
% r = 2;
% s = 2;
simVersion = "linear";
modes = 1;
pendTime = 0;
printDes = false;

% DEBUG
Td1
Td2
zeta1
zeta2

output = simInputShape(b,k,sys,sysRigid,Td,Td2,zeta,zeta2,tdes,tdessim,xEnd,...
        vStart,vEnd,fa11,fa12,fa21,fa22,p,q,r,s,st,shift,forwardF,...
        simVersion,modes,pendTime,printDes);
    
pos_sim     = output(:,1);
theta_sim   = output(:,2);
vel_sim     = output(:,3);
omega_sim   = output(:,4);
tsim = st:st:st*length(pos_sim);
texp = st:st:st*length(pos);
                
% Plot trimmed simulated and experimental velocity on same figure
figure();
plot(texp,vel,'LineWidth',2)
set(gca,'FontSize',16)
hold on;
plot(tsim,vel_sim,'LineWidth',2,'Color','#D95319')
hold off;
title('Experimental and Simulated Cart Velocities','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Velocity (m/s)','FontSize',18)


%% Generate system and input shaping impulses - 4 impulses
%{
b = 14.4;
k = 543.21;

[sys,Td1,Td2,zeta1,zeta2,overdamped] = sysCreate(b,k,false);

% Generate impulse times
t1      = 0;                                                            
t21     = Td1/2;                                                       
t22     = Td2/2;                                                       
toffset = t21+t22;
tf      = tdes-toffset;                                                 
t       = 0:st:tf;
tc      = 0:st:tdes;                                                    
lentc   = length(tc);

% Generate first mode impulse amplitudes. a and b are multiplicative
% factors to reproduce human subject execution error.
fa11    = 1;
fa12    = 1;
K1      = exp(-zeta1*pi/sqrt(1-zeta1^2));     
A11     = fa11*1/(1+K1);                                                   
A21     = fa12*K1/(1+K1);                                                  

% Generate second mode impulse amplitudes. c and d are multiplicative
% factors to reproduce human subject execution error.
fa21    = 1;
fa22    = 1;
K2      = exp(-zeta2*pi/sqrt(1-zeta2^2));     
A12     = fa21*1/(1+K2);                        
A22     = fa22*K2/(1+K2);             

% Create vector for first mode impulses
t21     = round(t21,digits(st));                                       
imp1    = zeros(1,length(t));                                           
imp1(int16(1+t1/st))    = A11;                                          
imp1(int16(1+t21/st))   = A21;                                          

% Create vector for second mode impulses
t22     = round(t22,digits(st));                                        
imp2    = zeros(1,length(t));                                           
imp2(int16(1+t1/st))    = A12;                                          
imp2(int16(1+t22/st))   = A22;                                          

% Convolve two sets of impulses
pulses  = conv(imp1,imp2);

% Find indicies of all pulses
pulsetimes = find(pulses); 
pulsePlot = pulses(pulsetimes);
pulsetimes = st*pulsetimes;

% Plot input shaping impulses
figure();
stem(pulsetimes,pulsePlot,'Marker','none','LineWidth',4)
set(gca,'FontSize',14)
title('Input Shaping Impulses - Full Internal Model','FontSize',20)
xlabel('Time (s)','FontSize',18)
ylabel('Impulse Amplitude','FontSize',18)
xlim([0 1])
ylim([0 1])
%}

%% Helper Function Definitions

% Returns order of magnitude of decimal number less than one
function y = digits(z)
    z = abs(z);
    y = 0;
    while (floor(z)~=z)
        y = y+1;
        z = z*10;
    end
end

