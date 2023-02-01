%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakRatios
%
% Takes a velocity profile, identifies the amplitudes of its local maxima,
% and calculates the asymmetry of the peaks.
% (Replaces peakAsymmetry.m for naming clarity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

%%%% SETTINGS
%%%% Settings: simulation type, block numbers, and trial numbers to analyze
% Choose range of simulation types to analyze
simStart = 1;
simEnd = 11;

% Choose which experimental blocks to analyze (all acceptable blocks
% would be 1-22)
blockStart  = 1;
blockEnd    = 22;

% Choose range of experimental trials to analyze
trialStart    = 1;
trialEnd      = 50;

%%%% Settings: plotting
plot = false;               % Plot trials (no plotting if false)
plotSave = true;            % Save figures
plotPeaks = true;           % Plot peak markers
plotSim = true;             % Plot simulated trials
plotExp = true;             % Plot experimental trials
plotDes = false;            % Plot desired input (ZFT)
tExp = [];
velExp = [];
velDes = [];
velChangeTimesExp = [];
velChangeAmpsExp = [];
plotHistogram = false;           % Plot histogram of peak ratios by subject

plotRows = 5;
plotCols = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIALIZATION
% Initialize arrays to hold ratios by trial, and median ratio by subject
velRatio = nan(trialStart-trialEnd+1, blockStart-blockEnd+1, simStart-simEnd+1);
medianRatio = nan(simEnd,blockEnd/2);

% Label for simulation type
simLabels = {'IS', 'MM-NF', 'SM-NF', 'FM-NF', 'RB-NF', 'NI-NF', ...
    'MM-FF', 'SM-FF', 'FM-FF', 'RB-FF', 'NI-FF'};

% Create table to hold median peak ratio values for each subject within a
% simulation
medianBySubject = nan(simEnd,blockEnd/2);
% {'Simulation Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
%     'S7'; 'S8'; 'S9'; 'S10'; 'S11'};
meanBySubject = nan(simEnd,blockEnd/2);

statsByBlock = cell(4,23);
statsByBlock(:,1) = {'Subject'; 'Mean Peak Ratio'; 'Median Peak Ratio';...
    'Std Dev Peak Ratio'}; 
statsByBlock(1,2:end) = {'S1_B3', 'S1_B4', 'S2_B3', 'S2_B4', 'S3_B3',...
    'S3_B4', 'S4_B3', 'S4_B4', 'S5_B3', 'S5_B4', 'S6_B3', 'S6_B4',...
    'S7_B3', 'S7_B4', 'S8_B3', 'S8_B4', 'S9_B3', 'S9_B4', 'S10_B3',...
    'S10_B4', 'S11_B3', 'S11_B4'};

% Create table to hold overall mean, median, and standard deviation values
% of velocity peak ratios for all simulations analyzed
overallResults = {'Simulation Type';...
                'Median Peak Ratio';...
                'Mean Peak Ratio';...
                'Std Dev of Peak Ratio';...
                'Mean of Median Peak Ratios';...
                'Std Dev of Median Peak Ratios';...
                'Significant to 5%?';
                '5% p value';...
                'Significant to 2%?';...
                '2% p value';...
                'Significant to 1%?';...
                '1% p value'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN LOOP
            
% Loop through simulation types
for simNum = simStart:simEnd
    
    % Get name of folder holding best-fit simulations for selected
    % simulation type
    parentFolder = simulationFolder(simNum);
    if parentFolder == "/"
        continue
    end
    
    % Loop through blocks
    for blockNum = blockStart:blockEnd
        
        % Get information about experimental block
        [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);
        
        % Determine if current block is block 3 or block 4
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

        if plot
            % Initialize figure window and assign handle
            figure('Name','Cart velocity','units','normalized','outerposition',[0 0 1 1]);
            cv = get(gcf,'Number');
            figTitle = ['Best Fit Cart Velocity Profiles: Subject ', subjNum, ' Block ', block3or4 ];
            sgtitle(figTitle,'FontSize',16);
        end
        
        % Loop through trials, calculating asymmetry of each one and
        % storing values in array
        for trialNum = trialStart:1:trialEnd
            
            % Skip invalid trials
            if ismember(trialNum,invalidTrials)
                velRatio(trialNum,blockNum,simNum) = NaN;
                continue
            end
            
            % Load best-fit simulated trial
            simProfiles = readmatrix(strcat(simFolder,...
                "/best-fit profiles/", num2str(trialNum), ".xlsx"));
            tSim        = simProfiles(:,1);
            velSim      = simProfiles(:,4);
            accSim      = simProfiles(:,6);
            st          = tSim(2) - tSim(1); % extract step size from time vector (only works if steps are consistent!)
            
            % Find velocity local minima and maxima
            velChangeTimesSim = find(diff(sign(accSim)))+1;
            
            % Combine pairs of changes where acceleration = 0
            flagInd = [];
            for ind = 1:length(velChangeTimesSim)
                if (accSim(velChangeTimesSim(ind)) == 0) && (accSim(velChangeTimesSim(ind+1)-1) == 0)
                    velChangeTimesSim(ind) = idivide(velChangeTimesSim(ind)+velChangeTimesSim(ind+1),int16(2));
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
%             medianBySubject(1,simNum+1) = {simLabels(simNum)};
%             medianBySubject(int8(blockNum/2),simNum) = ...
%                 {median(velRatio(:,blockNum-1:blockNum,simNum),'all','omitnan')};
            medianBySubject(simNum,int8(blockNum/2)) = ...
                median(velRatio(:,blockNum-1:blockNum,simNum),'all','omitnan');
            meanBySubject(simNum,int8(blockNum/2)) = ...
                mean(velRatio(:,blockNum-1:blockNum,simNum),'all','omitnan');
        end
        
        % Save figure contianing plots for block
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
%     [h5,p5,ci,stats] = ttest(cell2mat(medianBySubject(simNum+1,:)),1,...
%         'Tail',tail)
    [h5,p5,ci,stats] = ttest(medianBySubject(simNum+1,:),1,'Tail',tail)
    
%     [h2,p2,ci,stats] = ttest(cell2mat(medianBySubject(simNum+1,:)),1,...
%         'Tail',tail,'Alpha',.02)
    [h2,p2,ci,stats] = ttest(medianBySubject(simNum+1,:),1,'Tail',tail,'Alpha',.02)

%     [h1,p1,ci,stats] = ttest(cell2mat(medianBySubject(simNum+1,:)),1,...
%         'Tail',tail,'Alpha',.01)
    [h1,p1,ci,stats] = ttest(medianBySubject(simNum+1,:),1,'Tail',tail,'Alpha',.01)

    % Calculate t-value manually for sanity check
%     manualT = (mean(cell2mat(medianBySubject(2:end,simNum+1))) - 1)/...
%         std(cell2mat(medianBySubject(2:end,simNum+1)),1,'all','omitnan')/sqrt(10);
    
    % Put total results for simulation type into table
    overallResults(:,simNum-simStart+2) = ...
        {simLabels(simNum);...
        median(velRatio(:,:,simNum),'all','omitnan');...
        mean(velRatio(:,:,simNum),'all','omitnan');...
        std(velRatio(:,:,simNum),1,'all','omitnan');...
%         mean(cell2mat(medianBySubject(simNum+1,:)));...
%         std(cell2mat(medianBySubject(simNum+1,:)),1,'all','omitnan');...
        mean(medianBySubject(simNum+1,:));...
        std(medianBySubject(simNum+1,:),1,'all','omitnan');...
        h5; p5; h2; p2; h1; p1};
    
end

save('best fit simulation/peakRatioArrays.mat','medianBySubject','meanBySubject')

% % Rigid-Body Simplification Specifically
% rigidBodyVelRatios = velRatio(:,:,19);
% 
% filename = "VelRatios.xlsx";
% foldername = "best fit simulation/19. Rigid Body Simplification, Feedforward, B=6-10, K=75-125, 10k eval, 500ms delay/_Peak Ratios/";
% writematrix(rigidBodyVelRatios, strcat(foldername, filename))

%%
% Create histogram of peak ratios by subject
if plotHistogram

    % Initialize figure to hold subplots
    figure('Units','normalized','Position',[0,0,1,0.5]);
    subplotFig = tiledlayout(2,6);
    
    % Plot properties
    labelFontSize = 28;
    axisFontSize = 20;
    titleFontSize = 32;
    subtitleFontSize = 20;

    for blockNum = blockStart:blockEnd
        if ~mod(blockNum,2)

            [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
                minVelOutliers, invalidTrials] = blockDictionary(blockNum);

            subjVelRatios = velRatio(:,blockNum-1:blockNum,simNum);
            subjVelRatios = reshape(subjVelRatios,[],1);

            nexttile(subplotFig);
            histogram(subjVelRatios,'FaceColor',plotColor,'FaceAlpha',1);
            hold on;
            set(gca,'FontSize',axisFontSize)
            title(subjNum,'FontSize',subtitleFontSize);

       end
    end

    title(subplotFig,'Peak Ratios: Rigid-Body Simplification with Feedforward Force','FontSize',titleFontSize)
    xlabel(subplotFig,'Peak1/Peak2','FontSize',labelFontSize)
    ylabel(subplotFig,'Count','FontSize',labelFontSize)
    
end