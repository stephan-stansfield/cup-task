%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakRatios
%
% Takes a velocity profile, identifies the amplitudes of its local maxima,
% and calculates the ratio of the peaks. Calculates statistics across
% simulation type and subject.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

%%%% SETTINGS
%%%% Settings: simulation type, block numbers, and trial numbers to analyze
% Choose range of simulation types to analyze
simStart = 2;
simEnd = 6;

% Choose which experimental blocks to analyze (all acceptable blocks
% would be 1-22)
blockStart  = 1;
blockEnd    = 22;

% Choose range of experimental trials to analyze
trialStart    = 1;
trialEnd      = 50;

%%%% Settings: plotting
plot = false; % Plot trials
plotSave = true; % Save figures
plotPeaks = true; % Plot peak markers
plotSim = true; % Plot simulated trials
plotExp = true; % Plot experimental trials
plotDes = false; % Plot desired input (Zero Force Trajectory)
tExp = [];
velExp = [];
velDes = [];
velChangeTimesExp = [];
velChangeAmpsExp = [];

plotRows = 5;
plotCols = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIALIZATION
% Initialize arrays to hold ratios by trial, and mean & median ratio by subject
velRatio = nan(trialEnd-trialStart+1, blockEnd-blockStart+1, simEnd-simStart+2);
meanBySubject = nan(simEnd-simStart+1, int8(blockEnd-blockStart+1)/2);
medianBySubject = nan(simEnd-simStart+1, int8(blockEnd-blockStart+1)/2);

% multi-dimensional array to hold peak array for each trial across all
% internal model types:
peakRatioByTrial = nan(50, 22, 5);

% Label for simulation type
simLabels = {'EXP'; 'IS'; 'MM'; 'SM'; 'FM'; 'RB'; 'NI'};

% Create table to hold overall mean, median, and standard deviation values
% of velocity peak ratios for all simulations analyzed. The first row
% contains column header labels.
overallResults = {'Simulation Type',...
                'Mean Peak Ratio',...
                'Std Dev of Peak Ratios',...
                'Median Peak Ratio',...
                'Mean of Median Peak Ratios',...
                'Std Dev of Median Peak Ratios',...
                'Significant to 5%?',...
                '5% p value',...
                'Significant to 2%?',...
                '2% p value',...
                'Significant to 1%?',...
                '1% p value'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cycle through blocks
for blockNum = blockStart:blockEnd

    % Get information about experimental block
    [subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
        minVelOutliers, invalidTrials] = blockDictionary(blockNum);
    trialStr = erase(trialStr, "_trial_");
    
    % Determine if current block is block 3 or block 4
    block3or4 = blockNumFunc(blockNum);

    % Load file holding min vel and duration data
    load(strcat('peak times/', subjNum, '_B', block3or4))
    
    % Calculate velocity peak ratios for block of trials
    velRatio(:, blockNum, 1) = velAmpArray(:,1)./velAmpArray(:,3);
    
    % Compute per-subject statistics for both blocks 3 & 4
    if block3or4 == '4'
        meanBySubject(1, int8(blockNum / 2)) = ...
            mean(velRatio(:, blockNum - 1:blockNum, 1),'all','omitnan');
        medianBySubject(1, int8(blockNum / 2)) = ...
            median(velRatio(:, blockNum - 1:blockNum, 1),'all','omitnan');
    end

    % Conduct t-test to determine if mean of peak ratio medians is
    % significantly different from 1 at 5%, 2%, & 1%  signifiance levels
    tail = 'both';
    [h5,p5,ci5,stats] = ttest(medianBySubject(1,:),1,'Tail',tail)
    [h2,p2,ci2,stats] = ttest(medianBySubject(1,:),1,'Tail',tail,'Alpha',.02)
    [h1,p1,ci1,stats] = ttest(medianBySubject(1,:),1,'Tail',tail,'Alpha',.01)
    
    % Put overall results for experimental data into table
    overallResults(2, :) = ...
        {simLabels(1),...
        mean(velRatio(:,:,1),'all','omitnan'),...
        std(velRatio(:,:,1),1,'all','omitnan'),...
        median(velRatio(:,:,1),'all','omitnan'),...
        mean(medianBySubject(1,:)),...
        std(medianBySubject(1,:),1,'all','omitnan'),...
        h5, p5, h2, p2, h1, p1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulated Data - Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
% Loop through simulation types
for simNum = simStart:simEnd
    
    % Get name of folder holding best-fit simulations for selected
    % simulation type
    parentFolder = simulationFolder(simNum);
    if parentFolder == "/"
        continue
    end

    simNum = simNum + 1;
    
    % Loop through blocks
    for blockNum = blockStart:blockEnd
        
        % Get information about experimental block
        [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);
        
        % Determine if current block is block 3 or block 4
        block3or4 = blockNumFunc(blockNum);
    
        % Folder holding best-fit trials from current block
        if(isfolder(strcat("best fit simulation/",parentFolder,"S10_B3")))
            % Folders named by subject number & block number
            simFolder = strcat("best fit simulation/",parentFolder,subjNum,"_B",block3or4);
        else
            % Folders named based on subject initials & trial date (old convention):
            simFolder = strcat("best fit simulation/",parentFolder,subjStr,trialDate,trialStr);
            simFolder = erase(simFolder, "_trial_");
        end

        if plot
            % Initialize figure window and assign handle
            figure('Name','Cart velocity','units','normalized',...
                'outerposition',[0 0 1 1]);
            cv = get(gcf,'Number');
            figTitle = ['Best Fit Cart Velocity Profiles: Subject ', subjNum, ' Block ', block3or4 ];
            sgtitle(figTitle,'FontSize',16);
        end
        
        % Loop through trials, calculating asymmetry of each one and
        % storing values in array
        for trialNum = trialStart:1:trialEnd
            
            % Skip invalid trials
            if ismember(trialNum, invalidTrials)
                velRatio(trialNum, blockNum, simNum) = NaN;
                continue
            end
            
            % Load best-fit simulated trial
            simProfiles = readmatrix(strcat(simFolder,...
                "/best-fit profiles/", num2str(trialNum), ".xlsx"));
            tSim        = simProfiles(:,1);
            velSim      = simProfiles(:,4);
            accSim      = simProfiles(:,6);
            st          = tSim(2) - tSim(1);
            
            % Find velocity local minima and maxima
            velChangeTimesSim = find(diff(sign(accSim)))+1;
            
            % Combine pairs of changes where acceleration = 0
            flagInd = [];
            for ind = 1:length(velChangeTimesSim)
                if (accSim(velChangeTimesSim(ind)) == 0) &&...
                        (accSim(velChangeTimesSim(ind+1)-1) == 0)
                    velChangeTimesSim(ind) = idivide(velChangeTimesSim(ind)+...
                        velChangeTimesSim(ind+1),int16(2));
                    flagInd = [flagInd,ind];
                end
            end
            velChangeTimesSim(flagInd+1) = [];
            
            % Count number of velocity changes
            numPeaks = length(velChangeTimesSim);
            
            % Delete transient changes in vel from list of peaks & valleys
            velChangeTimesSim = peakFilter(numPeaks, velChangeTimesSim, velSim);
            
            % Velocity values at the local minima & maxima
            velChangeAmpsSim = velSim(velChangeTimesSim);
            
            % Calculate ratio of peak amplitudes if two peaks exist and
            % trial is not invalid
            if length(velChangeTimesSim) > 2
                peak1 = velSim(velChangeTimesSim(1));
                peak2 = velSim(velChangeTimesSim(3));
                velRatio(trialNum,blockNum,simNum) = peak1/peak2;
            else
                velRatio(trialNum,blockNum,simNum) = NaN;
            end
            
            if ismember(trialNum,invalidTrials)
                velRatio(trialNum,blockNum,simNum) = NaN;
                velChangeAmpsSim = [];
            end
            
            % Send trial and peak data to plotting function
            if plot
                figPlotAllTrials(cv, plotRows, plotCols, trialNum, trialStart,...
                   tExp, velExp, plotSim, plotExp, plotDes, tSim, velSim,...
                   plotPeaks, velChangeTimesSim, velChangeTimesExp, st,...
                   velChangeAmpsSim, velChangeAmpsExp, velDes) 
            end
        end
        % ^ End of looping through trials
        
        % Compute per-block statistics
        velRatioBlock = velRatio(:,blockNum,simNum);
        statsByBlock(2:4,blockNum+1) = {mean(velRatioBlock,'all','omitnan');...
            median(velRatioBlock,'all','omitnan');...
            std(velRatioBlock,1,'all','omitnan')};
        
        % Compute per-subject statistics for both blocks 3 & 4
        if ~mod(blockNum,2)
            medianBySubject(simNum, int8(blockNum / 2)) = ...
                median(velRatio(:, blockNum - 1:blockNum, simNum),'all','omitnan');
            meanBySubject(simNum, int8(blockNum / 2)) = ...
                mean(velRatio(:, blockNum - 1:blockNum, simNum),'all','omitnan');
        end
        
        % Save figure containing plots for block
        if plot && plotSave
            % Make output folder if it doesn't exist yet
            if(not(isfolder(strcat("best fit simulation/",parentFolder,"_Peak Ratios"))))
                mkdir(strcat("best fit simulation/",parentFolder),"_Peak Ratios");
            end
            
            saveFolder = strcat("best fit simulation/",parentFolder,"_Peak Ratios/");
            saveas(cv, strcat(saveFolder, "Sim Peaks_", subjNum, "_B", block3or4, ".png"));
        end
    end
    
    % Conduct t-test to determine if mean of peak ratio medians is
    % significantly different from 1 at 5%, 2%, & 1%  signifiance levels
    tail = 'both';
    [h5,p5,ci5,stats] = ttest(medianBySubject(simNum,:),1,'Tail',tail)
    [h2,p2,ci2,stats] = ttest(medianBySubject(simNum,:),1,'Tail',tail,'Alpha',.02)
    [h1,p1,ci1,stats] = ttest(medianBySubject(simNum,:),1,'Tail',tail,'Alpha',.01)
    
    % Put total results for simulation type into table
    overallResults(simNum-simStart+2, :) = ...
        {simLabels(simNum),...
        mean(velRatio(:,:,simNum),'all','omitnan'),...
        std(velRatio(:,:,simNum),1,'all','omitnan'),...
        median(velRatio(:,:,simNum),'all','omitnan'),...
        mean(medianBySubject(simNum,:)),...
        std(medianBySubject(simNum,:),1,'all','omitnan'),...
        h5, p5, h2, p2, h1, p1};
    
end

save('best fit simulation/peakRatioArrays.mat','medianBySubject','meanBySubject')