%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes data output from hapticmaster and plots variables against time.

clear all;
close all;
clc;

%{
% Initialize figure windows and assign variable names to figure numbers
figure('Name','Cart position','Units','inches','Position',[1,7.5,10,5.5])
cp = get(gcf,'Number');
title('Cart Position')
xlabel('time (s)')
ylabel('Cart Position (m)')

figure('Name','Ball angle','Units','inches','Position',[12.5,7.5,10,5.5])
bp = get(gcf,'Number');
title('Ball Angle')
xlabel('time (s)')
ylabel('Ball Angle (deg)')

figure('Name','Cart velocity','Units','inches','Position',[1,0,10,5.5])
cv = get(gcf,'Number');
title('Hand Velocity')
xlabel('time (s)')
ylabel('Hand Velocity (m/s)')

figure('Name','Ball angular velocity','Units','inches','Position',[12.5,0,10,5.5])
bv = get(gcf,'Number');
title('Ball Angular Velocity')
xlabel('time (s)')
ylabel('Ball Angular Velocity (deg/s)')

figure('Name','Cart acceleration','Units','inches','Position',[1,3.75,10,5.5]);
ca = get(gcf,'Number');
title('Cart Acceleration')
xlabel('time (s)')
ylabel('Cart Acceleration (m/s^2)')

figure('Name','Ball acceleration','Units','inches','Position',[12.5,3.75,10,5.5]);
ba = get(gcf,'Number');
title('Ball Acceleration')
xlabel('time (s)')
ylabel('Ball Acceleration (deg/s^2)')

figure('Name','Fx','Units','inches','Position',[1,1.75,10,5.5]);
fx = get(gcf,'Number')
title('Force, x-direction')
xlabel('time (s)')
ylabel('Force (N)')

figure('Name','Fy','Units','inches','Position',[12.5,1.75,10,5.5]);
fy = get(gcf,'Number')
title('Force, y-direction')
xlabel('time (s)')
ylabel('Force (N)')

figure('Name','Fz','Units','inches','Position',[1,0,10,5.5]);
fz = get(gcf,'Number')
title('Force, z-direction')
xlabel('time (s)')
ylabel('Force (N)')

figure('Name','Fball','Units','inches','Position',[12.5,0,10,5.5]);
fb = get(gcf,'Number')
title('Ball force')
xlabel('time (s)')
ylabel('Force (N)')
%}

% Plot position and velocity data for all trials within a range. Lines are
% added to indicate goal box position and time of entry.

% Load file with experimental data
load('LK_09Jul2019_11-02-22_trial_40.mat')

% Range of trial numbers to examine
numstart = 47;
numend = 49;

% Strings for specific subject, date, and trial
subjstr = 'LK';
datestr = '_09Jul2019';
trialstr = '_11-02-22_trial_';

for num = numstart:1:numend
    numstr = num2str(num);

    % Construct filename and load file of one experimental trial
    filestr = strcat(subjstr,datestr,trialstr,numstr);
    load(filestr)
    [pos,vel,acc,theta,omega,alpha,t] = trim(pos,vel,acc,theta,omega,alpha,t);
    
    % Find time when cart starts moving
    vmax            = max(vel);
    startIndex    = find(vel > 0.015*vmax);
    startTime       = t(startIndex(1));  
    
    % Find time when cart crosses into box
    enterIndex = find(pos > 0.23);
    if enterIndex
        enterTime = t(enterIndex(1));  
    else
        enterTime = 0;
    end
    
    % Find time when cart reaches center of box
    goalIndex = find(pos > 0.25);
    if goalIndex
        goalTime = t(goalIndex(1));  
    else
        goalTime = 0;
    end
    
    % Find time when cart crosses out of box
    exitIndex = find(pos > 0.27);
    if exitIndex
        exitTime = t(exitIndex(1));  
    else
        exitTime = 0;
    end
    
    figure();
    
    % Cart position
    subplot(2,2,1)
    plot(t,pos,'LineWidth',2)
    hold on
    line(t,0.23*ones(length(pos),1),'LineStyle','--','LineWidth',2,'Color','#EDB120')
    line(t,0.25*ones(length(pos),1),'LineStyle','--','LineWidth',2,'Color','#77AC30')
    line(t,0.27*ones(length(pos),1),'LineStyle','--','LineWidth',2,'Color','#A2142F')
    line(startTime*ones(100,1),linspace(0,0.3),'LineStyle','--','LineWidth',2,'Color','#7E2F8E')
    if enterTime
        line(enterTime*ones(100,1),linspace(0,0.3),'LineStyle','--','LineWidth',2,'Color','#EDB120')
    end
    if goalTime
        line(goalTime*ones(100,1),linspace(0,0.3),'LineStyle','--','LineWidth',2,'Color','#77AC30')
    end
    if exitTime
        line(exitTime*ones(100,1),linspace(0,0.3),'LineStyle','--','LineWidth',2,'Color','#A2142F')
    end
    title('Cart Position')
    hold off  

    % Ball angle
    subplot(2,2,2)
    plot(t,theta,'LineWidth',2)
    hold on
    line(t,zeros(length(pos),1),'LineStyle','--')
    line(startTime*ones(100,1),linspace(-10,10),'LineStyle','--','LineWidth',2,'Color','#7E2F8E')
    if enterTime
        line(enterTime*ones(100,1),linspace(-10,10),'LineStyle','--','LineWidth',2,'Color','#EDB120')
    end
    if goalTime
        line(goalTime*ones(100,1),linspace(-10,10),'LineStyle','--','LineWidth',2,'Color','#77AC30')
    end
    if exitTime
        line(exitTime*ones(100,1),linspace(-10,10),'LineStyle','--','LineWidth',2,'Color','#A2142F')
    end
    title('Ball Angle')
    hold off

    % Cart velocity
    subplot(2,2,3)
    plot(t,vel,'LineWidth',2)
    hold on
    line(startTime*ones(100,1),linspace(-0.1,0.5),'LineStyle','--','LineWidth',2,'Color','#7E2F8E')
    if enterTime
        line(enterTime*ones(100,1),linspace(-0.1,0.5),'LineStyle','--','LineWidth',2,'Color','#EDB120')
    end
    if goalTime
        line(goalTime*ones(100,1),linspace(-0.1,0.5),'LineStyle','--','LineWidth',2,'Color','#77AC30')
    end
    if exitTime
        line(exitTime*ones(100,1),linspace(-0.1,0.5),'LineStyle','--','LineWidth',2,'Color','#A2142F')
    end
    title('Cart Velocity')
    hold off

    % Ball angular velocity
    subplot(2,2,4)
    plot(t,omega,'LineWidth',2)
    hold on
    line(startTime*ones(100,1),linspace(-60,60),'LineStyle','--','LineWidth',2,'Color','#7E2F8E')
    if enterTime
        line(enterTime*ones(100,1),linspace(-60,60),'LineStyle','--','LineWidth',2,'Color','#EDB120')
    end
    if goalTime
        line(goalTime*ones(100,1),linspace(-60,60),'LineStyle','--','LineWidth',2,'Color','#77AC30')
    end
    if exitTime
        line(exitTime*ones(100,1),linspace(-60,60),'LineStyle','--','LineWidth',2,'Color','#A2142F')
    end
    title('Ball Velocity')
    hold off
    
    sgtitle(['Trial: ',num2str(num)],'FontSize',16)

end
    
%{
%%
% Force data
figure(fx)
plot(t,Fx)

figure(fy)
plot(t,Fy)

figure(fz)
plot(t,Fz)

figure(fb)
plot(t,Fball)

%%
%%%%% Change these variables for specific block and range of trials %%%%%


% Intialize variables to hold max parameter values
posmax      = 0;
thetamax    = 0;
velmax      = 0;
omegamax    = 0;
accmax      = 0;
alphamax    = 0;

% First trial in range to plot (trials start at 0)
numstart = 0;

% Last trial in range to plot (trials end at 49)
numend = 49;

% Number of rows and columns in subplot
plotrows = 5;
plotcols = 10;

% Strings for specific subject, date, and trial
subjstr = 'LK';
datestr = '_09Jul2019';
trialstr = '_11-02-22_trial_';

for num = numstart:1:numend
    numstr = num2str(num);

    % Construct filename and load file of one experimental trial
    filestr = strcat(subjstr,datestr,trialstr,numstr);
    load(filestr)
    
    %%%%% Treat data %%%%%
    % Offset cart position so it starts at x=0
    startpos = mean(pos(93:185));           % Average from 1s to 2s
    pos = pos - startpos;
    
    % Cut first 2 seconds of data
    t       = t(185:end);
    pos     = pos(185:end);
    theta   = theta(185:end);
    vel     = vel(185:end);
    omega   = omega(185:end);
    acc     = acc(185:end);
    
    % Find first and last points where velocity reaches 5% of maximum
    vmax        = max(vel);
    above5      = find(vel > 0.05*vmax);
    startindex  = above5(1);  
    stopindex   = above5(end);

    % Cut data before and after min velocity thresholds
    t       = t(startindex:stopindex);
    pos     = pos(startindex:stopindex);
    theta   = theta(startindex:stopindex);
    vel     = vel(startindex:stopindex);
    omega   = omega(startindex:stopindex);
    acc     = acc(startindex:stopindex);
    
    % Start time at 0 seconds
    t = t - t(1);
    
    % Convert angular quantities from radians to degrees
    theta = rad2deg(theta);
    omega = rad2deg(omega);  
    %%%%%
    
    % Plot current trial in subplots for each variable %
    % Cart position
    figure(cp);
    subplot(plotrows,plotcols,num-numstart+1)
    plot(t,pos,'LineWidth',2)
    title(['Trial' ' ' numstr]);
    xlabel('t')
    ylabel('cart pos')
    hold on;
    
    % Ball angle
    figure(bp);
    subplot(plotrows,plotcols,num-numstart+1)
    plot(t,theta,'LineWidth',2)
    title(['Trial' ' ' numstr]);
    xlabel('t')
    ylabel('angle')
    hold on;
    
    % Cart velocity
    figure(cv);
    subplot(plotrows,plotcols,num-numstart+1)
    plot(t,vel,'LineWidth',2)
    title(['Trial' ' ' numstr]);
    xlabel('t')
    ylabel('cart vel')
    hold on;
    
    % Ball angular velocity
    figure(bv);
    subplot(plotrows,plotcols,num-numstart+1)
    plot(t,omega,'LineWidth',2)
    title(['Trial' ' ' numstr]);
    xlabel('t')
    ylabel('ball vel')
    hold on;
    
    % Cart acceleration
    figure(ca);
    subplot(plotrows,plotcols,num-numstart+1)
    plot(t,acc,'LineWidth',2)
    title(['Trial' ' ' numstr]);
    xlabel('t')
    ylabel('cart acc')
    hold on;
    
    % If variable max for given trial is greater than previous trials,
    % assign to global max for that variable
    pmax = max(pos);
    tmax = max(theta);
    vmax = max(vel);
    omax = max(omega);
    amax = max(acc);
    
    if pmax > posmax
        posmax = pmax;
    end
    
    if tmax > thetamax
        thetamax = tmax;
    end
    
    if vmax > velmax
        velmax = vmax;
    end
    
    if omax > omegamax
        omegamax = omax;
    end
    
    if amax > accmax
        accmax = amax;
    end

end
    
% Adjust 'subjstr' for giving titles to figure windows
subjstr = strcat(subjstr,'\');

% Assign titles to figure windows
figure(cp);
figtitle = ['Cart Position Profiles,' ' ' subjstr datestr];
sgtitle(figtitle,'FontSize',16)
hold off;

figure(bp);
figtitle = ['Ball Angle Profiles,' ' ' subjstr datestr];
sgtitle(figtitle,'FontSize',16)
hold off;

figure(cv);
figtitle = ['Cart Velocity Profiles,' ' ' subjstr datestr];
sgtitle(figtitle,'FontSize',16)
hold off;

figure(bv);
figtitle = ['Ball Angular Velocity Profiles,' ' ' subjstr datestr];
sgtitle(figtitle,'FontSize',16)
hold off;

figure(ca);
figtitle = ['Cart Acceleration Profiles,' ' ' subjstr datestr];
sgtitle(figtitle,'FontSize',16)
hold off;

% Display global variable maxima for all trials
posmax
thetamax
velmax
omegamax
accmax
%}


function [pos,vel,acc,theta,omega,alpha,t] = trim(pos,vel,acc,theta,omega,alpha,t)

    % Offset position so it starts at x=0 and trim static section of data
    startpos    = mean(pos(1:93));                                          % Average position in first second
    pos         = pos - startpos;
    above3      = find(t > 3);
    startTime   = above3(1);  

    % Convert ball data to degrees from radians
    theta       = theta*360/(2*pi);
    omega       = omega*360/(2*pi);
    alpha       = alpha*360/(2*pi);

    pos         = pos(startTime:end);
    vel         = vel(startTime:end);
    acc         = acc(startTime:end);
    theta       = theta(startTime:end);
    omega       = omega(startTime:end);
    alpha       = alpha(startTime:end);
    t           = t(startTime:end);

end

