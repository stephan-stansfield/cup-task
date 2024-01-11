% figImpedanceScatter
% 
% Plots the mean best-fit impedance value by subject for a given simulation
% type in a histogram color-coded by subject.

clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data')); % DELETE?
saveFiles = false;

% Choose range of simulation types to analyze
simStart = 4;
simEnd = 4;

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
transparency = 0.65;

colorList = {[0 0 0], [0 0 0], [0 0 0], [0 0 0], [0 0 0]};

markerList = [  "diamond", ...  % Multi-mode input shaping
                'v', ...        % Slow mode internal model
                '^', ...        % Fast mode internal model
                "square", ...   % Rigid body internal model (hand impedance mode)
                'o'];           % No impedance internal model (pendulum mode)

simLabels = {'Multi-Mode', 'Slow Mode', 'Fast Mode', 'Rigid Body', 'No Impedance'};
xVals = categorical(simLabels,simLabels);

histogramK = figure();
histogramB = figure();

% Create color scheme
map = 'jet';
cmap = colormap(map);
rowmap = 1:23:256;

for simNum = simStart:simEnd
    
    % Folder within "best fit simulations" for this Control Model
    parentFolder = simulationFolder(simNum+1);  % note +1
    resultsTable = readtable(strcat("best fit simulation/",parentFolder,"Best-Fit Parameters.xlsx"));

    for subjNum = 1:11
        rangeStart = (subjNum-1)*100+1;
        rangeEnd = subjNum*100;
        meanArray(simNum, subjNum, 1) = mean(resultsTable{rangeStart:rangeEnd,'K'});
        meanArray(simNum, subjNum, 2) = mean(resultsTable{rangeStart:rangeEnd,'B'});
        medianArray(simNum, subjNum, 1) = median(resultsTable{rangeStart:rangeEnd,'K'});
        medianArray(simNum, subjNum, 2) = median(resultsTable{rangeStart:rangeEnd,'B'});

        % histograms
        blockNum = subjNum * 2;
        [~, ~, ~, ~, ~, plotColor, ~, ~] = blockDictionary(blockNum);
        subjK = resultsTable{rangeStart:rangeEnd,'K'};
        subjB = resultsTable{rangeStart:rangeEnd,'B'};

        set(groot,'CurrentFigure',histogramK);
        hold on;
        histogram(subjK, 'FaceColor', plotColor, 'EdgeColor', 'none', ...
            'FaceAlpha', transparency);
%         histogram(subjK, 'FaceColor', cmap(rowmap(subjNum), :),...
%             'EdgeColor', 'none', 'FaceAlpha', transparency);

        set(groot,'CurrentFigure',histogramB);
        hold on;
        histogram(subjB, 'FaceColor', plotColor, 'EdgeColor', 'none', ...
            'FaceAlpha', transparency);
%         histogram(subjB, 'FaceColor', cmap(rowmap(subjNum), :), ...
%             'EdgeColor', 'none', 'FaceAlpha', transparency);
    end

    if saveFiles
        mkdir(strcat("best fit simulation/",parentFolder),...
            "_Impedance Scatter Plots");
        saveFolder = strcat('best fit simulation/',parentFolder,...
            '_Impedance Scatter Plots'); 
    end
    
end

% Figure post-processing
% Means:
set(groot,'CurrentFigure',histogramK);
set(gca,'FontSize',axisFontSize)
xlabel('Best-Fit Stiffness by Subject (N/m)','FontSize',labelFontSize)
xlim([0 1000])
% legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', ...
%     'Location', 'northeast')
 
set(groot,'CurrentFigure',histogramB);
set(gca,'FontSize',axisFontSize)
xlabel('Best-Fit Damping by Subject (N-s/m)','FontSize',labelFontSize)
xlim([0 100])
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', ...
    'Location', 'northeast')
