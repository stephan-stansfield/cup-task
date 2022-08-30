%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trim and plot a single experimental trial that has been loaded into the
% workspace manually.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENTER BLOCK NUMBER FROM 1-22 AND TRIAL NUMBER FROM 1-50 HERE
blockNum = 22;
trialNum = 27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add relevant folders to the MATLAB path. Note that genpath adds the
% folder input and all of its subfolders
addpath('data', 'trim times', genpath('experimental data'));

% Get trial info
[subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
        minVelOutliers, invalidTrials] = blockDictionary(blockNum);
    
% Load start and stop indices of trials for this block
load(blockStr,'Expression1')

% Load data for this experimental trial
numStr = num2str(trialNum-1);
fileStr = strcat(subjStr,trialDate,trialStr,numStr);
load(fileStr,'pos','theta','vel','omega','acc','alpha','t')

% Load start and stop indices of this trial
start   = Expression1(trialNum,1);
stop    = Expression1(trialNum,2);

% Trim experimental data
st = 1;
[~,~,~,~,~,~,~,~,tExp,velRaw,accRaw] = trimData(pos,theta,vel,omega,...
    acc,alpha,t,st,start,stop);

% Block 3 or 4? (For figure title)
block3or4 = blockNumFunc(blockNum);

% Plot trimmed trial
figure();
plot(tExp,velRaw,"LineWidth",2)
title(strcat(subjNum, " Block ", block3or4, " Trial ", num2str(trialNum), " Trimmed Hand Velocity Profile"))
ylabel("Velocity (m/s)")
xlabel("Trimmed Time (s)")