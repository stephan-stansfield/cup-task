% minVelPlot
%
% Takes a block of kinematic trajectories and plots the values of the
% inter-peak velocity "dip" against the duration of the trajectory.
%
% This version of minVelPlot takes data on velocity minima and trimmed
% durations from saved .mat files.
%
% [Outdated: Simulated trajectories use the trimmed duration of the corresponding
% experimental trajectory for the duration value.]

clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data'));
addpath('peak times')

% Choose whether to plot simulated (true) or experimental (false) results
plotSim = true;

% Create tables to hold r and p values for each subject within a dataset
rBySubject = {'Data Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
    'S7'; 'S8'; 'S9'; 'S10'; 'S11'};
rBySubject(1,2) = {'Experimental'};
pBySubject = {'Data Type'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6';...
    'S7'; 'S8'; 'S9'; 'S10'; 'S11'};
pBySubject(1,2) = {'Experimental'};

% Choose whether to plot by both blocks for each subject (true) or by
% each block individually (false)
pairs = true;

% Choose whether to save figures (Note: "pairs" must also be true)
save = true;

if save
    saveFolder = strcat("best fit simulation/00. Experimental Trials/_Experimental Min Velocity Correlation/");
    pngSuffix = ".png";
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

    % Load file holding min vel and duration data
    load(strcat('peak times/',subjNum,'_B',block3or4))
    
    % Remove trials that don't have values for minimum velocity (invalid
    % trials and trials without two clear peaks)
    durationArray = rmmissing(durationArray);
    velAmpArray = rmmissing(velAmpArray);
    
    % Grab velocity minima (second column of array)
    velMinBlock = velAmpArray(:,2);
    
    % Add loaded values to array holding both blocks for one subject
    if block3or4 == '3'
        correlationArray = [durationArray, velMinBlock];
    else
        correlationArray = [correlationArray;
                            durationArray, velMinBlock];
    end

    % If grouping both blocks and on block 3, move on to block 4
    if pairs && block3or4 == '3'
        continue
    end

    % Fit linear trendline to data points
    durations = correlationArray(:,1);
    minVels = correlationArray(:,2);
    [fitPoly,S] = polyfit(durations, minVels, 1);

%         % Calculate fitted line R^2 value
%         Rsq = 1 - (S.normr/norm(minVels - mean(minVels)))^2;

    % Calculate Pearson's correlation coefficient, r, and significance
    % level, p.
    [R,P] = corrcoef(durations,minVels);
    r = R(1,2);
    p = P(1,2);

    % Print block information and linear fit values
    if pairs
        blockStr = subjNum;
    end
    disp(['Block: ', blockStr]);
    disp(['Line of best fit: y = ', num2str(fitPoly(1)), '*x + ', num2str(fitPoly(2))]);
    disp(['r value: ', num2str(r)]);
    disp(['p value: ', num2str(p)]);

    rBySubject(int8(blockNum/2)+1,2) = {r};
    pBySubject(int8(blockNum/2)+1,2) = {p};

    % Find the max and min durations, round them up and down to the
    % nearest 0.1 s, respectively
    maxDuration = ceil(max(durations)*10)/10;
    minDuration = floor(min(durations)*10)/10;

    % For plotting, create array with fitted values within plot region
    xfit = linspace(minDuration,maxDuration);
    yfit = fitPoly(1)*xfit + fitPoly(2);

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

    if pairs && save
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

title(subplotFig,'Human Subject Data','FontSize',titleFontSize)
xlabel(subplotFig,'Movement Duration (s)','FontSize',labelFontSize);
ylabel(subplotFig,'Inte-peak Minimum Velocity (m/s)','FontSize',labelFontSize);


if pairs && save
    saveas(gcf, strcat(saveFolder, "subplots.png"))          
end

close all;


