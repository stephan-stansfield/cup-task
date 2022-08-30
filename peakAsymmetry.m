%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peakAsymmetry
%
% Takes a velocity profile, identifies the amplitudes of its local maxima,
% and calculates the asymmetry of the peaks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTINGS
%%%% Settings: simulation type, block numbers, and trial numbers to analyze
% Choose range of simulation types to analyze
simStart = 25;
simEnd = 25;

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
% Initialize array to hold velocity ratios
velRatio = nan(trialStart-trialEnd+1, blockStart-blockEnd+1, simStart-simEnd+1);

% Label for simulation type
simLabels = {'IS', 'MM-FF', 'NI-FF', 'RB-FF', 'SM-FF', 'FM-FF', 'MM-NF',...
    'NI-NF', 'RB-NF', 'SM-NF', 'FM-NF', 'SB', 'SB-FF', 'SB-NF', 'RB-FF-FI',...
    'RB-FF-NI-10K', 'RB-FF-K-10K', 'RB-FF-B6.10-K75.125-10K',...
    'RB-FF-B6.10-K75.125-10K-dur', 'RB-FF-B6.10-K75.175', 'delay250-625',...
    'FBO','dur_corr','RB-FF-B6.20-K0.250-10K','RB-FF-B7.20-K0.250-10K',...
    '26. IS delay500', '27. MM delay500','28. RB-FF-B7.20-K0.250-10K-delay250',...
    '29. RB-FF-B0.50-K0.750-10K-delay500'};

% Create table to hold median peak ratio values for each subject within a
% simulation
medianBySubject = {'Simulation Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
    'S7'; 'S8'; 'S9'; 'S10'; 'S11'};

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
                '2% p value'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN LOOP
            
% Cycle through simulation types
for simNum = simStart:simEnd
    
    % 1. Nominal Input Shaping (No Impedance Model)
    % 2. Multi-Mode Input Shaping, FF
    % 3. No Impedance IS, FF
    % 4. Rigid Body IS, FF
    % 5. Slow Mode IS, FF
    % 6. Fast Mode IS, FF
    % 7. Multi-Mode Input Shaping, No FF
    % 8. No Impedance IS, No FF
    % 9. Rigid Body IS, No FF
    % 10. Slow Mode IS, No FF
    % 11. Fast Mode IS, No FF
    % 12. Submovements, No Impedance
    % 13. Submovements, Impedance, FF
    % 14. Submovements, Impedance, No FF
    
    % Get name of folder holding best-fit simulations for selected
    % simulation type
    parentFolder = simulationFolder(simNum);
    if parentFolder == "/"
        continue
    end
    
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
            
%             % DEBUG
%             disp('Trial:')
%             trialNum
            
%             % Skip invalid trials and vel vs. duration correlation outliers
%             if ismember(trialNum,invalidTrials) || ismember(trialNum,minVelOutliers)
%                 velRatio(trialNum,blockNum,simNum) = NaN;
%                 continue
%             end
            
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
            medianBySubject(1,simNum+1) = {simLabels(simNum)};
            medianBySubject(int8(blockNum/2)+1,simNum+1) = ...
                {median(velRatio(:,blockNum-1:blockNum,simNum),'all','omitnan')};
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
    % significantly different from 1 at 5% and 2% signifiance levels
    tail = 'both';
    [h5,p5,ci,stats] = ttest(cell2mat(medianBySubject(2:12,simNum+1)),1,...
        'Tail',tail)
    
    [h2,p2,ci,stats] = ttest(cell2mat(medianBySubject(2:12,simNum+1)),1,...
        'Tail',tail,'Alpha',.02)

    % Calculate t-value manually for sanity check
%     manualT = (mean(cell2mat(medianBySubject(2:end,simNum+1))) - 1)/...
%         std(cell2mat(medianBySubject(2:end,simNum+1)),1,'all','omitnan')/sqrt(10);
    
    % Put total results for simulation type into table
    overallResults(:,simNum-simStart+2) = {simLabels(simNum);...
        median(velRatio(:,:,simNum),'all','omitnan');...
        mean(velRatio(:,:,simNum),'all','omitnan');...
        std(velRatio(:,:,simNum),1,'all','omitnan');...
        mean(cell2mat(medianBySubject(2:end,simNum+1)));...
        std(cell2mat(medianBySubject(2:end,simNum+1)),1,'all','omitnan');...
        h5; p5; h2; p2};
    
end

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