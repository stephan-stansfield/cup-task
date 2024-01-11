% figImpedanceScatter
% 
% Plots the mean best-fit impedance value by subject for different
% simulation types

clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data')); % DELETE?
saveFiles = false;

% Choose range of simulation types to analyze
simStart = 1;
simEnd = 5;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Array to hold results. Dimension lengths are # of simulations analyzed, #
% of subjects (11), and # of parameters (2: K & B)
meanArray = zeros(simStart-simEnd+1, 11, 2);
medianArray = zeros(simStart-simEnd+1, 11, 2);

% Plotting settings
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;
marSz = 150;                        % Size of scatter plot markers
blue = [0, 0.4470, 0.7410];         % MATLAB default blue
green = [0.4660 0.6740 0.1880];     % MATLAB default green
transparency = 0.75;

colorList = {[0 0 0], [0 0 0], [0 0 0], [0 0 0], [0 0 0]};

markerList = [  "diamond", ...  % Multi-mode input shaping
                'v', ...        % Slow mode internal model
                '^', ...        % Fast mode internal model
                "square", ...   % Rigid body internal model (hand impedance mode)
                'o'];           % No impedance internal model (pendulum mode)

simLabels = {'Multi-Mode', 'Slow Mode', 'Fast Mode', 'Rigid Body', 'No Impedance'};
xVals = categorical(simLabels,simLabels);

meanScatterK = figure();
meanScatterB = figure();
% medianScatterK = figure();
% medianScatterB = figure();

for simNum = simStart:simEnd
    
    % Folder within "best fit simulations" for this Control Model
    parentFolder = simulationFolder(simNum+1);  % note +1
    resultsTable = readtable(strcat("best fit simulation/",parentFolder,"Best-Fit Parameters.xlsx"));

    for subjNum = 1:11
        [~, ~, ~, ~, ~, plotColor, ~, ~] = blockDictionary(2*subjNum);
        colorList{subjNum} = plotColor;

        rangeStart = (subjNum-1)*100+1;
        rangeEnd = subjNum*100;
        meanArray(simNum, subjNum, 1) = mean(resultsTable{rangeStart:rangeEnd,'K'});
        meanArray(simNum, subjNum, 2) = mean(resultsTable{rangeStart:rangeEnd,'B'});
        medianArray(simNum, subjNum, 1) = median(resultsTable{rangeStart:rangeEnd,'K'});
        medianArray(simNum, subjNum, 2) = median(resultsTable{rangeStart:rangeEnd,'B'});

        % Plot means
        set(groot,'CurrentFigure',meanScatterK);
        Kscatter = scatter(xVals(simNum), meanArray(simNum,subjNum,1),marSz,'Marker',markerList(simNum),...
            'MarkerEdgeColor','none', 'MarkerFaceColor',plotColor);
        hold on;
        alpha(transparency)

        set(groot,'CurrentFigure',meanScatterB);
        Bscatter = scatter(xVals(simNum), meanArray(simNum,subjNum,2),marSz,'Marker',markerList(simNum),...
            'MarkerEdgeColor','none', 'MarkerFaceColor',plotColor);
        hold on;
        alpha(transparency)
        
    end

    if saveFiles
        mkdir(strcat("best fit simulation/",parentFolder),...
            "_Impedance Scatter Plots");
        saveFolder = strcat('best fit simulation/',parentFolder,...
            '_Impedance Scatter Plots'); 
    end
    
%     % Means:
%     set(groot,'CurrentFigure',meanScatterK);
%     Kscatter = scatter(xVals(simNum), meanArray(simNum,:,1),marSz,'Marker',markerList(simNum),...
%         'MarkerEdgeColor','none', 'MarkerFaceColor',colorList{simNum});
%     hold on;
%     alpha(transparency)
% %     Kscatter.MarkerFaceAlpha = transparency*ones(11,1);
% %     Kscatter.MarkerEdgeAlpha = transparency*ones(1,11);
%     
% 
%     set(groot,'CurrentFigure',meanScatterB);
%     Bscatter = scatter(xVals(simNum), meanArray(simNum,:,2),marSz,'Marker',markerList(simNum),...
%         'MarkerEdgeColor','none', 'MarkerFaceColor',colorList{simNum});
%     alpha(transparency)
%     hold on;

    % Medians:
%     set(groot,'CurrentFigure',medianScatterK);
%     scatter(xVals(simNum), medianArray(simNum,:,1),marSz,'Marker',markerList(simNum),...
%         'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
%     hold on;
% 
%     set(groot,'CurrentFigure',medianScatterB);
%     scatter(xVals(simNum), medianArray(simNum,:,2),marSz,'Marker',markerList(simNum),...
%         'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
%     hold on;
    
end

% Figure post-processing
% Means:
set(groot,'CurrentFigure',meanScatterK);
set(gca,'FontSize',axisFontSize)
% title('Per-Subject Mean Best-Fit Stiffness by Control Model','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
ylabel('Best-Fit Stiffness (N/m)','FontSize',labelFontSize)
ylim([0, 1000])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
 
set(groot,'CurrentFigure',meanScatterB);
set(gca,'FontSize',axisFontSize)
% title('Per-Subject Mean Best-Fit Damping by Control Model','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
ylabel('Best-Fit Damping (N-s/m)','FontSize',labelFontSize)
ylim([0, 100])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
% legend([markerObjects{1} markerObjects{2} markerObjects{3} markerObjects{4},...
%     markerObjects{5}], xVals(1:5),'Location','northeast')


% % Medians:
% set(groot,'CurrentFigure',medianScatterK);
% title('Per-Subject Median Best-Fit Stiffness by Control Model','FontSize',titleFontSize,'FontWeight','Normal')
% hold on;
% set(gca,'FontSize',axisFontSize)
% ylabel('Stiffness (N/m)','FontSize',labelFontSize)
% ylim([0, 800])
% set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
%  
% set(groot,'CurrentFigure',medianScatterB);
% title('Per-Subject Median Best-Fit Damping by Control Model','FontSize',titleFontSize,'FontWeight','Normal')
% hold on;
% set(gca,'FontSize',axisFontSize)
% ylabel('Damping (N-s/m)','FontSize',labelFontSize)
% ylim([5, 50])
% set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
