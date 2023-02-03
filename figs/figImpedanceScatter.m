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
simStart = 1; % note off by one to leave out input shaping w/ no impedance
simEnd = 10;

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

colorList = {   blue, blue, blue, blue, blue,... % No feedforward
                green, green, green, green, green};     % Feedforward

markerList = [  "diamond",...                           % Multi-mode input shaping
                'v',...                                 % Slow mode internal model
                '^',...                                 % Fast mode internal model
                "square",...                            % Rigid body internal model (hand impedance mode)
                'o',...                                 % No impedance internal model (pendulum mode)
                "diamond", 'v', '^', "square", 'o'];    % Repeat for feedforward

simLabels = {'MM-NF', 'SM-NF', 'FM-NF', 'RB-NF', 'NI-NF', ...
             'MM-FF', 'SM-FF', 'FM-FF', 'RB-FF', 'NI-FF'};
xVals = categorical(simLabels,simLabels);

meanScatterK = figure();
meanScatterB = figure();
medianScatterK = figure();
medianScatterB = figure();

for simNum = simStart:simEnd
    
    % Folder within "best fit simulations" for this simulation type
    parentFolder = simulationFolder(simNum+1);  % note +1
    resultsTable = readtable(strcat("best fit simulation/",parentFolder,"Best-Fit Parameters.xlsx"));

    for subjNum = 1:11
        rangeStart = (subjNum-1)*100+1;
        rangeEnd = subjNum*100;
        meanArray(simNum, subjNum, 1) = mean(resultsTable{rangeStart:rangeEnd,'K'});
        meanArray(simNum, subjNum, 2) = mean(resultsTable{rangeStart:rangeEnd,'B'});
        medianArray(simNum, subjNum, 1) = median(resultsTable{rangeStart:rangeEnd,'K'});
        medianArray(simNum, subjNum, 2) = median(resultsTable{rangeStart:rangeEnd,'B'});
    end

    if saveFiles
        mkdir(strcat("best fit simulation/",parentFolder),...
            "_Impedance Scatter Plots");
        saveFolder = strcat('best fit simulation/',parentFolder,...
            '_Impedance Scatter Plots'); 
    end
    
    set(groot,'CurrentFigure',meanScatterK);
    scatter(xVals(simNum), meanArray(simNum,:,1),marSz,'Marker',markerList(simNum),...
        'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
    hold on;

    set(groot,'CurrentFigure',meanScatterB);
    scatter(xVals(simNum), meanArray(simNum,:,2),marSz,'Marker',markerList(simNum),...
        'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
    hold on;

    set(groot,'CurrentFigure',medianScatterK);
    scatter(xVals(simNum), medianArray(simNum,:,1),marSz,'Marker',markerList(simNum),...
        'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
    hold on;

    set(groot,'CurrentFigure',medianScatterB);
    scatter(xVals(simNum), medianArray(simNum,:,2),marSz,'Marker',markerList(simNum),...
        'MarkerEdgeColor',colorList{simNum}, 'MarkerFaceColor',colorList{simNum});
    hold on;
    
end

% Figure post-processing
set(groot,'CurrentFigure',meanScatterK);
title('Per-Subject Mean Best-Fit Stiffness by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
set(gca,'FontSize',axisFontSize)
ylabel('Stiffness (N/m)','FontSize',labelFontSize)
ylim([0, 800])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
 
set(groot,'CurrentFigure',meanScatterB);
title('Per-Subject Mean Best-Fit Damping by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
set(gca,'FontSize',axisFontSize)
ylabel('Damping (N-s/m)','FontSize',labelFontSize)
ylim([5, 50])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead

set(groot,'CurrentFigure',medianScatterK);
title('Per-Subject Median Best-Fit Stiffness by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
set(gca,'FontSize',axisFontSize)
ylabel('Stiffness (N/m)','FontSize',labelFontSize)
ylim([0, 800])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
 
set(groot,'CurrentFigure',medianScatterB);
title('Per-Subject Median Best-Fit Damping by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal')
hold on;
set(gca,'FontSize',axisFontSize)
ylabel('Damping (N-s/m)','FontSize',labelFontSize)
ylim([5, 50])
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead
