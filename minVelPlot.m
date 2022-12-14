% minVelPlot
%
% Takes a block of kinematic trajectories and plots the values of the
% inter-peak velocity "dip" against the duration of the trajectory.

% [Outdated: Simulated trajectories use the trimmed duration of the corresponding
% experimental trajectory for the duration value.]

clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data'));

% Choose whether to plot simulated (true) or experimental (false) results
plotSim = true;

% Label for data type (experimental or type of simulation)
dataLabels = {'IS', 'MM-FF', 'NI-FF', 'RB-FF', 'SM-FF', 'FM-FF', 'MM-NF',...
    'NI-NF', 'RB-NF', 'SM-NF', 'FM-NF', 'SB', 'SB-FF', 'SB-NF', 'RB-FF-FI', ...
    'RB-FF-NI-10K', 'RB-FF-K-10K', 'RB-FF-B6.10-K75.125-10K', 'Test', ...
    'RB-FF-B6.10-K75.125-10K_corr', 'RB-FF-B6.10-K75.175', 'delay250-625',...
    'FBO', 'dur_corr','RB-FF-B6.20-K0.250-10K', 'RB-FF-B7.20-K0.250-10K',...
    '26. IS delay500', '27. MM delay500','28. RB-FF-B7.20-K0.250-10K-delay250',...
    '29. RB-FF-B0.50-K0.750-10K-delay500', 'Exp'};

% Create tables to hold r and p values for each subject within a dataset
rBySubject = {'Data Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
    'S7'; 'S8'; 'S9'; 'S10'; 'S11'};
pBySubject = {'Data Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
    'S7'; 'S8'; 'S9'; 'S10'; 'S11'};

% Choose range of simulation types to analyze
simStart = 28;
simEnd = 29;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Choose whether to save figures (Note: "pairs" must also be true)
saveFiles = true;

% Choose whether to plot by both blocks for each subject (true) or by
% each block individually (false)
pairs = true;

% Time step of simulated profile. NOTE that this has to match the step
% used when simulating the profile. Usually this is 1 ms.
st = 0.001;

% Name of folder within "best fit simulations" folder holding blocks of 
% interest
for simNum = simStart:simEnd
    
    parentFolder = simulationFolder(simNum);

    if saveFiles
        if plotSim
            mkdir(strcat("best fit simulation/",parentFolder),...
                "_Minimum Velocity Plots");
            saveFolder = strcat('best fit simulation/',parentFolder,...
                '_Minimum Velocity Plots/'); 
        else
            saveFolder = strcat('best fit simulation/',...
                '_Experimental Min Velocity/');
            pngSuffix = '.png';
        end
        
        if simNum == 99
            saveFolder = strcat('best fit simulation/',...
                '_Experimental Min Velocity/');
            pngSuffix = '.png';
        end
            
    end
    
    % Initialize figure to hold subplots
    figure('Units','normalized','Position',[0,0,1,0.5]);
    subplotFig = tiledlayout(2,6);


     % Choose which experimental blocks to visualize
    for blockNum = 1:22

        % Get information about current block.
        [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);
        
        block3or4 = blockNumFunc(blockNum);
        
        % Folder holding best-fit trials from current block
        if(isfolder(strcat("best fit simulation/",parentFolder,"S1_B3")))
            % Folders named by subject number & block number
            simFolder = strcat("best fit simulation/",parentFolder,subjNum,"_B",block3or4);
        else
            % Folders named based on subject initials & trial date (old convention):
            simFolder = strcat("best fit simulation/",parentFolder,subjStr,trialDate,trialStr);
            simFolder = erase(simFolder, "_trial_");
        end

        % Load start and stop indices of experimental trials for this block
        load(blockStr,'Expression1')

        if pairs
            if block3or4 == '3'
                % Initialize array to hold minimum velocity and trial duration values
                minVelArray = nan(100,2);

                % Initialize array to hold time of minimum velocity values
                timeArray = nan(100,1);
            end
        else
            % Initialize array to hold minimum velocity and trial duration values
            minVelArray = nan(50,2);

            % Initialize array to hold time of minimum velocity values
            timeArray = nan(50,1);
        end

        % Loop through trials, store min velocity and duration values in array
        for num = numStart:1:numEnd

            % DEBUG
%             disp('Trial #:')
%             num
            
            % Skip invalid trials and trials that did not have two clear peaks
            if ismember(num,minVelOutliers)
                % DEBUG
%                 disp('Outlier!')
                continue
            end

            numStr = num2str(num-1);

            % Construct file name and load file of one experimental trial
            fileStr = strcat(subjStr,trialDate,trialStr,numStr);
            load(fileStr,'pos','theta','vel','omega','acc','alpha','t')

            % Load start & stop indices of corresponding experimental trial
            start   = Expression1(num,1);
            stop    = Expression1(num,2);

            % Trim experimental data to start & stop indices
            [pos,theta,vel,omega,acc,alpha,tdes] = trimData(pos,theta,vel,...
                omega,acc,alpha,t,st,start,stop);

            % Load best-fit SIMULATED trial
            simProfiles = readmatrix(strcat(simFolder, "/best-fit profiles/", num2str(num), ".xlsx"));
            tdesSim     = simProfiles(end,1);
            posSim      = simProfiles(:,2);
            thetaSim    = simProfiles(:,3);
            velSim      = simProfiles(:,4);
            omegaSim    = simProfiles(:,5);
            accSim      = simProfiles(:,6);
            alphaSim    = simProfiles(:,7);
            dur_corr    = simProfiles(1,8);

            if plotSim
                % Find minimum velocity of simulated profile
                [timeDip, valDip] = minVelFind(velSim,accSim,st);
            else
                % Find minimum velocity of experimental profile
                [timeDip, valDip] = minVelFind(vel,acc,st);
            end
            
            if simNum == 99
                % Find minimum velocity of experimental profile
                [timeDip, valDip] = minVelFind(vel,acc,st);
            end

%             % Skip values that were flagged (valDip = NaN)
%             if isnan(timeDip)
%                 continue
%             end

            % Put minimum velocity and experimental duration into array
            if pairs
                if block3or4 == '4'
                    arrayNum = num+50;
                else
                    arrayNum = num;
                end
            else
                arrayNum = num;
            end

            % Superseded: use experimental trial duration for calculating
            % simulated duration-min. vel correlation
%             % Put corresponding durations and velocity dip values in array
%             minVelArray(arrayNum,:) = [tdes, valDip];
            
            % Put corresponding durations and velocity dip values in array
            minVelArray(arrayNum,:) = [dur_corr, valDip];

            % Put time of minimum velocity into separate array to check work
            timeArray(arrayNum) = timeDip;

        end

        % If grouping both blocks and on block 3, move on to block 4
        if pairs && (block3or4 == '3')
            continue
        end

        % Delete flagged elements from array
%         [zeroRow, ~] = find(~minVelArray);
%         minVelArray(zeroRow,:) = [];
        [flagRow, ~] = find(isnan(minVelArray));
        minVelArray(flagRow,:) = [];

        % Fit linear trendline to data points
        durations = minVelArray(:,1);
        minVels = minVelArray(:,2);
        [poly,S] = polyfit(durations, minVels, 1);

%         % Calculate fitted line R^2 value
%         Rsq = 1 - (S.normr/norm(minVels - mean(minVels)))^2;

        % Calculate Pearson's correlation coefficient, r, and significance
        % level, p
        [R,P] = corrcoef(durations,minVels);
        r = R(1,2);
        p = P(1,2);

        % Print block information and linear fit values
        if pairs
            blockStr = subjNum;
        end
        disp(['Block: ', blockStr]);
        disp(['Line of best fit: y = ', num2str(poly(1)), '*x + ', num2str(poly(2))]);
        disp(['r value: ', num2str(r)]);
        disp(['p value: ', num2str(p)]);

        rBySubject(1,simNum+1) = {dataLabels(simNum)};
        rBySubject(int8(blockNum/2)+1,simNum+1) = {r};
        
        pBySubject(1,simNum+1) = {dataLabels(simNum)};
        pBySubject(int8(blockNum/2)+1,simNum+1) = {p};

        % Find the max and min durations, round them up and down to the
        % nearest 0.1 s, respectively
        maxDuration = ceil(max(durations)*10)/10;
        minDuration = floor(min(durations)*10)/10;

        % Create array with linear fit values within plot region
        xfit = linspace(minDuration,maxDuration);
        yfit = poly(1)*xfit + poly(2);

        % Figure properties
        labelFontSize = 28;
        axisFontSize = 20;
        titleFontSize = 32;
        subtitleFontSize = 20;
%         labelFontSize = 16;
%         axisFontSize = 16;
%         titleFontSize = 16;
        figTitle = blockStr;

        % Plot value of minimum velocity versus trial duration in
        % individual figures
        figure();
        plot(durations,minVels,'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
        hold on;
        plot(xfit,yfit,'LineWidth',2,'Color',plotColor)
        set(gca,'FontSize',axisFontSize)
        xlabel('Trial Duration (s)','FontSize',labelFontSize)
        ylabel('Minimum Velocity (m/s)','FontSize',labelFontSize)
        title(figTitle,'FontSize',titleFontSize);
        xlim([minDuration, maxDuration]);
        hold off;

        if pairs && saveFiles
            saveas(gcf, strcat(saveFolder, subjNum, ".png"))          
        end
        
        % Plot correlations for all subjects in one figure
        nexttile(subplotFig);
        if blockNum/2 == 1
            nexttile(1)
        end
        plot(durations,minVels,'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
        hold on;
        plot(xfit,yfit,'LineWidth',2,'Color',plotColor)
        set(gca,'FontSize',axisFontSize)
        title(figTitle,'FontSize',subtitleFontSize);
        xlim([minDuration, maxDuration]);
        
    end

    
    title(subplotFig,{'Inter-Peak Minimum Velocity vs. Duration','Rigid-Body Simplification with Feedforward Force'},'FontSize',titleFontSize)
    xlabel(subplotFig,'Movement Duration (s)','FontSize',labelFontSize);
    ylabel(subplotFig,'Inter-Peak Minimum Velocity (m/s)','FontSize',labelFontSize);
    
    if pairs && saveFiles
        saveas(gcf, strcat(saveFolder, "min vel vs duration correlation.png"))          
        save(strcat(saveFolder, 'minVelData.mat'), 'durations', 'minVels')
    end

end


