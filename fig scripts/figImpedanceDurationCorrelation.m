% figImpedanceDurationCorrelation
% 
% Plots the mean best-fit impedance value versus trial duration by subject 
% for different simulation types

clear all;
close all;
clc;

% addpath('data', 'trim times', genpath('experimental data')); % DELETE?
saveFiles = true;

% Choose range of simulation types to analyze
simStart = 5;
simEnd = 5;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Arrays to hold results. Dimension lengths are # of simulations analyzed, 
% # of subjects (11), and # of parameters (2: K & B)
bestFitStiffArray = zeros(simStart-simEnd+1, 11, 100);
bestFitDurArray = zeros(simStart-simEnd+1, 11, 100);

% Plotting settings
titleFontSize = 26;
subtitleFontSize = 22;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;
marSz = 150;                        % Size of scatter plot markers
black = [0 0 0];
orange = [0.9290 0.6940 0.1250];    % MATLAB default orange
green = [0.4660 0.6740 0.1880];     % MATLAB default green

colorList = {   black, orange, orange, orange, orange, orange,...   % No feedforward
                green, green, green, green, green};                 % Feedforward

markerList = [  "pentagram",...                         % Original input shaping
                "diamond",...                           % Multi-mode input shaping without feedforward
                "v",...                                 % Slow mode internal model w/out ff
                "^",...                                 % Fast mode internal model w/out ff
                "square",...                            % Rigid body internal model (hand impedance mode) w/out ff
                "o",...                                 % No impedance internal model (pendulum mode) w/out ff
                "diamond", "v", "^", "square","o"];     % Repeat internal models with feedforward

rsquareArray = cell(11, 11);
sseArray = cell(11, 11);

for simNum = simStart:simEnd

    % Initialize figure to hold subplots
    figure('Units','normalized','Position',[0,0,1,0.5]);
    subplotFig = tiledlayout(2,6);
    
    % Folder within "best fit simulations" for this simulation type
    parentFolder = simulationFolder(simNum);
    resultsTable = readtable(strcat("best fit simulation/",parentFolder,"Best-Fit Parameters.xlsx"));

    if saveFiles
        mkdir(strcat("best fit simulation/",parentFolder),...
            "_Impedance-Duration Correlation Plots");
        saveFolder = strcat('best fit simulation/',parentFolder,...
            '_Impedance-Duration Correlation Plots'); 
    end
    
    for subjNum = 1:11
        rangeStart = (subjNum-1)*100+1;
        rangeEnd = subjNum*100;
        bestFitStiffArray(simNum, subjNum, 1:100) = resultsTable{rangeStart:rangeEnd,'K'};
%         bestFitArray(simNum, subjNum, 2) = mean(resultsTable{rangeStart:rangeEnd,'B'});

        for trial = 1:50
            blockFolder = strcat("S", num2str(subjNum), "_B3/best-fit profiles/");
            trialTable = readtable(strcat("best fit simulation/",parentFolder,blockFolder,num2str(trial),".xlsx"));
            bestFitDurArray(simNum,subjNum,trial) = trialTable{end,'Time_s_'};
        end
        for trial = 51:100
            blockFolder = strcat("S", num2str(subjNum), "_B4/best-fit profiles/");
            trialTable = readtable(strcat("best fit simulation/",parentFolder,blockFolder,num2str(trial-50),".xlsx"));
            bestFitDurArray(simNum,subjNum,trial) = trialTable{end,'Time_s_'};
        end

        durVector = bestFitDurArray(simNum, subjNum, :);
        stiffVector = bestFitStiffArray(simNum, subjNum, :);
        durVector = reshape(durVector, [numel(durVector), 1]);
        stiffVector = reshape(stiffVector, [numel(stiffVector), 1]);

        excludedInd = (durVector > 2.5) + (stiffVector > 300) + (stiffVector < 5); % exclude outliers
        excludedInd(excludedInd > 1) = 1;
        [f, gof] = fit(durVector, stiffVector, 'exp1', 'Exclude', excludedInd);
        rsquareArray{simNum, subjNum} = gof.rsquare;
        sseArray{simNum, subjNum} = gof.sse;

        figure();
        [~, ~, ~, ~, ~, plotColor, ~, ~] = blockDictionary(subjNum*2);
%         s = scatter(reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
%             reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]),...
%             marSz,'Marker',markerList(simNum),'MarkerEdgeColor',plotColor,...
%             'MarkerFaceColor',plotColor);
%         hold on;
%         plot(f)
%         set(s, 'LineWidth', 2)
%         f = plot(f,reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
%             reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]));
        scatter(reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
            reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]),...
            marSz,'Marker',markerList(simNum),'MarkerEdgeColor',plotColor,...
            'MarkerFaceColor',plotColor);
        hold on;
%         plot(reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
%             reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]), ...
%             'LineStyle','none');
        fitFig = plot(f);
        set(fitFig, 'LineWidth', 3, 'LineStyle','--', 'Color','k')
        legend('off')
        titleName = {'Best-Fit Stiffness vs. Trial Duration',...
            strcat(parentFolder(5:end-1), " S", num2str(subjNum))};
        title(titleName,'FontSize',titleFontSize,'FontWeight','Normal')
        set(gca,'FontSize',axisFontSize,'Color','none')
        ylabel('Stiffness (N/m)','FontSize',labelFontSize)
        xlim([0.75 2.5])
        ylim([0, 300])
        xlabel('Trial Duration (s)','FontSize',labelFontSize)

        saveas(gcf, strcat(saveFolder, "/S", num2str(subjNum), ".png"))
%         close % DEBUG

        % Plot correlations for all subjects in one figure
        nexttile(subplotFig);
        scatter(reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
            reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]),...
            marSz,'Marker',markerList(simNum),'MarkerEdgeColor',plotColor,...
            'MarkerFaceColor',plotColor);
        hold on;
        fitFig = plot(f);
        set(fitFig, 'LineWidth', 3, 'LineStyle','--', 'Color','k')
        legend('off')
        set(gca,'FontSize',axisFontSize,'Color','none')
        xlabel([])
        ylabel([])
        xlim([0.75 2.5])
        ylim([0, 300])
        title(strcat("S", num2str(subjNum)),'FontSize',subtitleFontSize);

    end

%     title(subplotFig, {'Best-Fit Stiffness vs. Trial Duration',...
%             strcat(parentFolder(5:end-1))}, 'FontSize', titleFontSize);
    xlabel(subplotFig, 'Trial Duration (s)','FontSize',labelFontSize);
    ylabel(subplotFig, 'Stiffness (N/m)','FontSize',labelFontSize);
    saveas(gcf, strcat(saveFolder, "/impedance-duration correlation.png"))

    close
    
end

%%
% Sandbox to work on figure settings after data has been generated

% figure();
% [~, ~, ~, ~, ~, plotColor, ~, ~] = blockDictionary(subjNum*2);
% % myFig = plot(f, durVector, stiffVector, marSz,'Marker',markerList(simNum),'MarkerEdgeColor',plotColor,...
% %     'MarkerFaceColor',plotColor);
% myFig = plot(f, durVector, stiffVector);
% % set(myFig, 'lineWidth', 2)
% legend('off')
% hold on;
% % scatter(reshape(bestFitDurArray(simNum,subjNum,1:100), [1 100]), ...
% %     reshape(bestFitStiffArray(simNum,subjNum,1:100), [1 100]),...
% %     marSz,'Marker',markerList(simNum),'MarkerEdgeColor',plotColor,...
% %     'MarkerFaceColor',plotColor,'LineWidth',0.5);
% 
% titleName = {'Best-Fit Stiffness vs. Trial Duration',...
%     strcat(parentFolder(5:end-1), " S", num2str(subjNum))};
% title(titleName,'FontSize',titleFontSize,'FontWeight','Normal')
% set(gca,'FontSize',axisFontSize)
% ylabel('Stiffness (N/m)','FontSize',labelFontSize)
% xlim([0.75 2.5])
% ylim([0, 300])
% xlabel('Trial Duration (s)','FontSize',labelFontSize)
