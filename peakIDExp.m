%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakIDExp.m: Takes the experimental data, runs the peak filtering
% algorithm, and saves a file containing the identified local minima and
% maxima for each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%%%% SETTINGS
% Choose which experimental blocks to analyze (all acceptable blocks
% would be 1-22)
blockStart  = 1;
blockEnd    = 22;

% Choose range of experimental trials to analyze
trialStart    = 1;
trialEnd      = 50;

plotRows = 5;
plotCols = 10;

st = 0.001; % Not used, but vestigial input for trimData function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% INITIALIZATION
% Add relevant folders to the MATLAB path. Note that genpath adds the
% folder input and all of its subfolders
addpath('data', 'trim times', genpath('experimental data'));

% Initialize array to hold velocity ratios
velRatio = nan(trialEnd-trialStart+1, blockEnd-blockStart+1);

% Create table to hold median peak ratio values for each subject within a
% simulation
medianBySubject = {'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';'S7'; 'S8'; 'S9';...
    'S10'; 'S11'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN LOOP

% Cycle through blocks
for blockNum = blockStart:blockEnd

    % Get information about experimental block
    [subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
        minVelOutliers, invalidTrials] = blockDictionary(blockNum);
    
    % Load start and stop indices of trials for this block
    load(blockStr,'Expression1')
    
    % Create arrays for this block to hold saved data
    velAmpArray = NaN(50, 3);
    velTimeArray = NaN(50, 3);
    durationArray = NaN(50, 1);

    % Loop through trials, calculating hand velocity local minima and
    % maxima of each one and storing values in array
    for trialNum = trialStart:1:trialEnd
        
        numStr = num2str(trialNum-1);

        % Skip invalid trials and trials without two clear peaks
        if ismember(trialNum,invalidTrials) || ismember(trialNum,minVelOutliers)
            velRatio(trialNum,blockNum) = NaN;
            continue
        end

        % Construct file name and load file of one experimental trial
        fileStr = strcat(subjStr,trialDate,trialStr,numStr);
        load(fileStr,'pos','theta','vel','omega','acc','alpha','t')
        
        % Load start and stop indices of corresponding trial
        start   = Expression1(trialNum,1);
        stop    = Expression1(trialNum,2);

        % Trim experimental data
        [~,~,~,~,~,~,~,~,tExp,velRaw,accRaw] = trimData(pos,theta,vel,omega,...
            acc,alpha,t,st,start,stop);

        % Find velocity local minima and maxima
        velChangeIndicesExp = find(diff(sign(accRaw)))+1;

        % Combine pairs of changes where acceleration = 0
        flagInd = [];
        for ind = 1:length(velChangeIndicesExp)
            if (accRaw(velChangeIndicesExp(ind)) == 0) && (accRaw(velChangeIndicesExp(ind+1)-1) == 0)
                velChangeIndicesExp(ind) = idivide(velChangeIndicesExp(ind)+velChangeIndicesExp(ind+1),int16(2));
                flagInd = [flagInd,ind];
            end
        end
        velChangeIndicesExp(flagInd+1) = [];

        % Count number of velocity changes
        numPeaks = length(velChangeIndicesExp);

        % Delete transient changes in velocity from list of peaks & valleys
        velChangeIndicesExp = peakFilterExp(numPeaks, velChangeIndicesExp, velRaw, tExp);
        
        % Velocity values at the local minima & maxima
        velChangeAmpsExp = velRaw(velChangeIndicesExp);
        velChangeTimesExp = tExp(velChangeIndicesExp);
        
        if length(velChangeAmpsExp) > 3
            velChangeAmpsExp = velChangeAmpsExp(1:3);
            velChangeTimesExp = velChangeTimesExp(1:3);
        end
        
        % If less than 3 peaks, fill in rest of values with NaNs
        % Maybe not necessary
        
        % Trial duration
        trialDuration = t(start) - t(stop);
        
        % Save trial duration, velocity amplitude, and velocity timing
        % values to arrays
        velAmpArray(trialNum,:) = velChangeAmpsExp;
        velTimeArray(trialNum,:) = velChangeTimesExp;
        durationArray(trialNum) = t(stop) - t(start);

    end

    % Get current block number (3 or 4) as string
    block3or4 = blockNumFunc(blockNum);
    
    % Compute per-subject statistics for both blocks 3 & 4
    if block3or4 == '4'
        medianBySubject(int8(blockNum/2),2) = ...
            {median(velRatio(:,blockNum-1:blockNum),'all','omitnan')};
    end
    
    % Save peak, valley, and duration variables to .mat file
    saveFolder = strcat("best fit simulation/Experimental Trials/Peak and Duration Data/");
    saveFile = strcat(subjNum, "_B", block3or4, ".mat");
    fileName = strcat(saveFolder,saveFile);
    save(fileName, 'velAmpArray', 'velTimeArray', 'durationArray');

end