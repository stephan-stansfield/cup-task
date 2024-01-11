% Print experimental trajectories and different best-fit trajectories on
% top of each other


%% 4-Impulse, No Timing Variation
close all;
clear all;
clc;

folder = 'plot data/4 Impulse No Timing Var/';

load([folder,'C.mat'])
load([folder,'avgArray.mat'])
load([folder,'4ImpulseF.mat'],'vel_sim')

vel_4ImpF = vel_sim;

load([folder,'4ImpulseNoF.mat'],'vel_sim')

vel_4ImpNoF = vel_sim;

lenExp = length(vel);
len4ImpF = length(vel_4ImpF);
len4ImpNoF = length(vel_4ImpNoF);

maxLen = max([lenExp, len4ImpF, len4ImpNoF]);

% Pad arrays that are shorter with zeros at the end
if lenExp < maxLen
    diff = maxLen - lenExp;
    vel = padarray(vel, [0,diff],'post');
end

if len4ImpF < maxLen
    diff = maxLen - len4ImpF;
    vel_4ImpF = padarray(vel_4ImpF,[diff,0],'post');
end

if len4ImpNoF < maxLen
    diff = maxLen - len4ImpNoF;
    vel_4ImpNoF = padarray(vel_4ImpNoF,[diff,0],'post');
end

for num = 1:50
    array = C{num};
    trialLength = max(size(array));
    diff = maxLen - trialLength;
    
    if diff > 0
        paddedArray = zeros(6,maxLen);
        paddedArray(1,:) = padarray(array(1,:), [0,diff],'replicate','post');
        paddedArray(2,:) = padarray(array(2,:), [0,diff],'replicate','post');
        paddedArray(3:6,:) = padarray(array(3:6,:), [0,diff],'post');

       % Replace array with padded version in cell
       C{num} = paddedArray;
    end
end
    
% Define time vector
st = 0.001;
t = 0:st:(maxLen-1)*.001;

figure();
for i = 1:50                                                                % Plot individual trials
    array = C{i};
    plot(t,array(3,:),'Color','#d6d6d6')
    hold on;
end
plot(t,vel,'Color','#000000','LineWidth',2)                                 % Experimental Average
plot(t,vel_4ImpF,'Color','#A2142F','LineWidth',2)                           % With F
plot(t,vel_4ImpNoF,'Color','#4DBEEE','LineWidth',2)                         % Wihtout F
title('4-Impulse Input Shaping, Feedforward F vs. No Feedforward F')
xlabel('Time (s)')
ylabel('Hand Velocity (m/s)')


%% 4-Impulse, Timing Variation
close all;
clear all;
clc;

folder = 'plot data/4 Impulse Timing Var/';

load([folder,'C.mat'])
load([folder,'avgArray.mat'])
load([folder,'4ImpulseF.mat'],'vel_sim')

vel_4ImpF = vel_sim;

load([folder,'4ImpulseNoF.mat'],'vel_sim')

vel_4ImpNoF = vel_sim;

lenExp = length(vel);
len4ImpF = length(vel_4ImpF);
len4ImpNoF = length(vel_4ImpNoF);

maxLen = max([lenExp, len4ImpF, len4ImpNoF]);

% Pad arrays that are shorter with zeros at the end
if lenExp < maxLen
    diff = maxLen - lenExp;
    vel = padarray(vel, [0,diff],'post');
end

if len4ImpF < maxLen
    diff = maxLen - len4ImpF;
    vel_4ImpF = padarray(vel_4ImpF,[diff,0],'post');
end

if len4ImpNoF < maxLen
    diff = maxLen - len4ImpNoF;
    vel_4ImpNoF = padarray(vel_4ImpNoF,[diff,0],'post');
end

for num = 1:50
    array = C{num};
    trialLength = max(size(array));
    diff = maxLen - trialLength;
    
    if diff > 0
        paddedArray = zeros(6,maxLen);
        paddedArray(1,:) = padarray(array(1,:), [0,diff],'replicate','post');
        paddedArray(2,:) = padarray(array(2,:), [0,diff],'replicate','post');
        paddedArray(3:6,:) = padarray(array(3:6,:), [0,diff],'post');

       % Replace array with padded version in cell
       C{num} = paddedArray;
    end
end
    
% Define time vector
st = 0.001;
t = 0:st:(maxLen-1)*.001;

figure();

for i = 1:50                                                                % Plot individual trials
    array = C{i};
    plot(t,array(3,:),'Color','#d6d6d6','HandleVisibility','off')
    hold on;
end
plot(t,vel,'Color','#000000','LineWidth',2)                                 % Experimental Average
plot(t,vel_4ImpF,'Color','#A2142F','LineWidth',2)                           % With F
plot(t,vel_4ImpNoF,'Color','#4DBEEE','LineWidth',2)                         % Wihtout F
title('Feedforward F vs. No Feedforward F','FontSize',24)
xlabel('Time (s)','FontSize',22,'FontWeight','bold')
ylabel('Hand Velocity (m/s)','FontSize',22,'FontWeight','bold')
legend('Exp. Avg.','With F','With No F','Location','northeast')
ax = gca;
ax.FontSize = 20; 

%% 4-Impulse vs 2-Impulse, Timing Variation
close all;
clear all;
clc;

folder = 'plot data/4 vs 2 Impulse Timing Var/';

load([folder,'C.mat'])
load([folder,'avgArray.mat'])
load([folder,'4ImpulseF.mat'],'vel_sim')

vel_4ImpF = vel_sim;

load([folder,'2ImpulseF.mat'],'vel_sim')

vel_2ImpF = vel_sim;

lenExp = length(vel);
len4ImpF = length(vel_4ImpF);
len2ImpF = length(vel_2ImpF);

maxLen = max([lenExp, len4ImpF, len2ImpF]);

% Pad arrays that are shorter with zeros at the end
if lenExp < maxLen
    diff = maxLen - lenExp;
    vel = padarray(vel, [0,diff],'post');
end

if len4ImpF < maxLen
    diff = maxLen - len4ImpF;
    vel_4ImpF = padarray(vel_4ImpF,[diff,0],'post');
end

if len2ImpF < maxLen
    diff = maxLen - len2ImpF;
    vel_2ImpF = padarray(vel_2ImpF,[diff,0],'post');
end

for num = 1:50
    array = C{num};
    trialLength = max(size(array));
    diff = maxLen - trialLength;
    
    if diff > 0
        paddedArray = zeros(6,maxLen);
        paddedArray(1,:) = padarray(array(1,:), [0,diff],'replicate','post');
        paddedArray(2,:) = padarray(array(2,:), [0,diff],'replicate','post');
        paddedArray(3:6,:) = padarray(array(3:6,:), [0,diff],'post');

       % Replace array with padded version in cell
       C{num} = paddedArray;
    end
end
    
% Define time vector
st = 0.001;
t = 0:st:(maxLen-1)*.001;

figure();

for i = 1:50                                                                % Plot individual trials
    array = C{i};
    plot(t,array(3,:),'Color','#d6d6d6','HandleVisibility','off')
    hold on;
end
plot(t,vel,'Color','#000000','LineWidth',2)                                 % Experimental Average
plot(t,vel_4ImpF,'Color','#A2142F','LineWidth',2)                           % 4 impulse
plot(t,vel_2ImpF,'Color','#4DBEEE','LineWidth',2)                           % 2 impulse
title('4- vs. 2-Impulse Input Shaping','FontSize',24)
xlabel('Time (s)','FontSize',22,'FontWeight','bold')
ylabel('Hand Velocity (m/s)','FontSize',22,'FontWeight','bold')
ylim([0,0.35])
legend('Exp. Avg.','4-Imp.','2-Imp.','Location','northeast')
ax = gca;
ax.FontSize = 20; 

%% Timing Variation vs. No Timing Variation, 4-Impulse
close all;
clear all;
clc;

folder = 'plot data/Timing Var vs No Timing Var, 4 Impulse/';

load([folder,'C.mat'])
load([folder,'avgArray.mat'])
load([folder,'4ImpulseF.mat'],'vel_sim')

vel_4ImpF = vel_sim;

load([folder,'4ImpulseFNoTime.mat'],'vel_sim')

vel_4ImpFNoTime = vel_sim;

lenExp = length(vel);
len4ImpF = length(vel_4ImpF);
len4ImpFNoTime = length(vel_4ImpFNoTime);

maxLen = max([lenExp, len4ImpF, len4ImpFNoTime]);

% Pad arrays that are shorter with zeros at the end
if lenExp < maxLen
    diff = maxLen - lenExp;
    vel = padarray(vel, [0,diff],'post');
end

if len4ImpF < maxLen
    diff = maxLen - len4ImpF;
    vel_4ImpF = padarray(vel_4ImpF,[diff,0],'post');
end

if len4ImpFNoTime < maxLen
    diff = maxLen - len4ImpFNoTime;
    vel_4ImpFNoTime = padarray(vel_4ImpFNoTime,[diff,0],'post');
end

for num = 1:50
    array = C{num};
    trialLength = max(size(array));
    diff = maxLen - trialLength;
    
    if diff > 0
        paddedArray = zeros(6,maxLen);
        paddedArray(1,:) = padarray(array(1,:), [0,diff],'replicate','post');
        paddedArray(2,:) = padarray(array(2,:), [0,diff],'replicate','post');
        paddedArray(3:6,:) = padarray(array(3:6,:), [0,diff],'post');

       % Replace array with padded version in cell
       C{num} = paddedArray;
    end
end
    
% Define time vector
st = 0.001;
t = 0:st:(maxLen-1)*.001;

figure();

for i = 1:50                                                                % Plot individual trials
    array = C{i};
    plot(t,array(3,:),'Color','#d6d6d6','HandleVisibility','off')
    hold on;
end
plot(t,vel,'Color','#000000','LineWidth',2)                                 % Experimental Average
plot(t,vel_4ImpF,'Color','#A2142F','LineWidth',2)                           % Timing Var
plot(t,vel_4ImpFNoTime,'Color','#4DBEEE','LineWidth',2)                     % No Timing Var

title('Impulse Timing Variability vs. No Variability','FontSize',24)
xlabel('Time (s)','FontSize',22,'FontWeight','bold')
ylabel('Hand Velocity (m/s)','FontSize',22,'FontWeight','bold')
legend('Exp. Avg.','Timing Var','No Timing Var','Location','northeast')
ax = gca;
ax.FontSize = 20; 
