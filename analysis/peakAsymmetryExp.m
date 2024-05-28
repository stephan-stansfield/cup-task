%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakAsymmetryExp
%
% Takes a velocity profile, identifies the amplitudes of its local maxima,
% and calculates the asymmetry of the peaks. This script is modified to
% analyze the experimentally gathered data, not simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%%%% SETTINGS
% Choose which experimental blocks to analyze (all subjects blocks 3 & 4
% would be 1-22)
blockStart = 1;
blockEnd = 22;

% Choose range of experimental trials to analyze
trialStart = 1;
trialEnd = 50;

%%%% Settings: plotting
plotTrial = false; % Plot trials (no plotting if false)
plotSave = false; % Save figures
plotPeaks = false; % Plot peak markers
plotSim = false; % Plot simulated trials
plotExp = false; % Plot experimental trials
plotDes = false; % Plot desired input (Zero Force Trajectory)
tSim = [];
velSim = [];
velDes = [];
velChangeTimesSim = [];
velChangeAmpsSim = [];

plotRows = 5;
plotCols = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% INITIALIZATION
% Add relevant folders to the MATLAB path
addpath('peak times');

% Initialize array to hold velocity ratios
velRatio = nan(trialEnd-trialStart+1, blockEnd-blockStart+1);

% Create table to hold median peak ratio values for each subject within a
% simulation
meanBySubject = {'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';'S7'; 'S8'; 'S9';...
    'S10'; 'S11'};
meanBySubject = nan(1,blockEnd/2);
medianBySubject = {'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';'S7'; 'S8'; 'S9';...
    'S10'; 'S11'};
medianBySubject = nan(1,blockEnd/2);

% Create table to hold overall mean, median, and standard deviation values
% of velocity peak ratios for all simulations analyzed
overallResults = {'Data Type';...
                'Median Peak Ratio';...
                'Mean Peak Ratio';...
                'Std Dev of Peak Ratio';...
                'Mean of Median Peak Ratios';...
                'Std Dev of Median Peak Ratios';...
                'Significant to 5%?';
                '5% p value';...
                'Significant to 1%?';...
                '1% p value'};
            
st = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN LOOP

% Cycle through blocks
for blockNum = blockStart:blockEnd

    % Get information about experimental block
    [subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
        minVelOutliers, invalidTrials] = blockDictionary(blockNum);
    trialStr = erase(trialStr, "_trial_");
    
    % Determine if current block is block 3 or block 4
    block3or4 = blockNumFunc(blockNum);
    
    if plotTrial    
        % Initialize figure window and assign handle
        figure('Name','Hand velocity','units','normalized','outerposition',...
            [0 0 1 1]);
        cv = get(gcf,'Number');
        figTitle = ['Experimental Hand Velocity Profiles: Subject ',subjNum,...
            ' Block ',block3or4 ];
        sgtitle(figTitle,'FontSize',16);
    end

    % Load file holding min vel and duration data
    load(strcat('peak times/',subjNum,'_B',block3or4))
    
    velRatio(:,blockNum) = velAmpArray(:,1)./velAmpArray(:,3);
    
    if plotTrial
        % Loop through trials to plot them with peak markers
        for trialNum = trialStart:1:trialEnd
            % Load trimmed experimental trial
            fileStr = strcat("experimental data trimmed/",subjStr,...
                trialDate,trialStr,"/",num2str(trialNum), ".mat");
            load(fileStr,'pos','theta','vel','omega','acc','alpha','tdes')
            t = 0:st:tdes;
            velChangeTimesExp = velTimeArray(trialNum,:);
            velChangeAmpsExp = velAmpArray(trialNum,:);
    
            % Send trial and peak data to plotting function
            figPlotAllTrials(cv, plotRows, plotCols, trialNum, trialStart,...
               t, vel, plotSim, plotExp, plotDes, tSim, velSim,...
               plotPeaks, velChangeTimesSim, velChangeTimesExp, st,...
               velChangeAmpsSim, velChangeAmpsExp, velDes) 
        end
    end
    
    % Compute per-subject statistics for both blocks 3 & 4
    if block3or4 == '4'
        meanBySubject(int8(blockNum/2)) = ...
            mean(velRatio(:,blockNum-1:blockNum),'all','omitnan');
        medianBySubject(int8(blockNum/2)) = ...
            median(velRatio(:,blockNum-1:blockNum),'all','omitnan');
    end
    
    % Save figure contianing plots for block
    if plotTrial && plotSave
        saveFolder = strcat("best fit simulation/0. Experimental Trials/_Peak Markers/");
        saveas(cv, strcat(saveFolder, "Sim Peaks_", subjNum, "_B", block3or4, ".png"));
    end

end

% Calculate standard deviation of subject means
stdMeanBySubject = std(meanBySubject);
