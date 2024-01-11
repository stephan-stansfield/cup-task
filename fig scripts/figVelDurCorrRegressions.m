% figVelocityDurationCorrelation (previously minVelPlot)
%
% Takes a block of kinematic trajectories and plots the linear regression 
% of the values of the interpeak velocity "dip" against the duration of 
% the trajectory for each subject.

clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data'));

% Choose whether to plot simulated (true) or experimental (false) results
plots = true; % whether to plot figures at all
plotSim = true;

% Choose whether to save data & figures (Note: "pairs" must also be true)
saveFiles = true;

% Choose whether to plot by both blocks for each subject (true) or by
% each block individually (false)
pairs = true;

% Choose range of simulation types to analyze
simStart = 1;
simEnd = 6;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Titles for subplots
subplotTitle = {'Nominal Input Shaping', 'Multi-Mode', 'Slow Mode', ...
    'Fast Mode', 'Rigid Body', 'No Impedance'};

% Initialize figure to hold subplots
if plots
    figure('Units','normalized','Position',[0,0,0.5,1]);
    subplotFig = tiledlayout(3,2);
end

for simNum = simStart:simEnd
    
    nexttile(subplotFig);

    % Folder within "best fit simulations" for this simulation type
    parentFolder = simulationFolder(simNum);

    % Load r and p values
        load(strcat('best fit simulation/', parentFolder, ...
            '_Minimum Velocity Plots/minVelData.mat'));

    if saveFiles
        if plotSim
            mkdir(strcat("best fit simulation/",parentFolder), ...
                "_Minimum Velocity Plots");
            saveFolder = strcat('best fit simulation/',parentFolder, ...
                '_Minimum Velocity Plots/'); 
        else
            saveFolder = strcat('best fit simulation/', ...
                '_Experimental Min Velocity/');
            pngSuffix = '.png';
        end
        
        if simNum == 99
            saveFolder = strcat('best fit simulation/', ...
                '_Experimental Min Velocity/');
            pngSuffix = '.png';
        end
            
    end

     % Choose which experimental blocks to visualize
    for blockNum = 1:22

        block3or4 = blockNumFunc(blockNum);

        % If grouping both blocks and on block 3, move on to block 4
        if pairs && (block3or4 == '3')
            continue
        end

        % If on block 4, plot regression line

        % Get information about current block
        [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);

        % Get r and p for current subject
        r = rArray{int8(blockNum/2)};
        p = pArray{int8(blockNum/2)};
        poly = polyArray{int8(blockNum/2)};

        % Print block information and linear fit values
        if pairs
            blockStr = subjNum;
        end
        disp(['Block: ', blockStr]);
        disp(['Line of best fit: y = ', num2str(poly(1)), '*x + ', num2str(poly(2))]);
        disp(['r value: ', num2str(r)]);
        disp(['p value: ', num2str(p)]);

        % Find the max and min durations, round them up and down to the
        % nearest 0.1 s, respectively
        maxDuration = 1.6;
        minDuration = 1.1;

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

        if plots            
            % Plot correlations for all subjects in one figure
            if blockNum/2 == 1
                continue
            end
%             set(groot, 'CurrentFigure', subplotFig)
            plot(xfit,yfit,'LineWidth',2,'Color',plotColor)
            hold on;
            set(gca,'FontSize',axisFontSize)
            title(subplotTitle{simNum},'FontSize',subtitleFontSize);
%             th = title(figTitle,'FontSize',subtitleFontSize);
%             titlePos = get(th, 'position')
%             titlePos(2) = titlePos(2) - 1
%             set(th, 'position', titlePos);
%             xlim([minDuration, maxDuration]);
        end
        
    end

    if plots
%         set(groot, 'CurrentFigure', subplotFig)
        xlabel(subplotFig, 'Movement Duration (s)','FontSize',labelFontSize);
        ylabel(subplotFig, 'Interpeak Minimum Velocity (m/s)','FontSize',labelFontSize);
        
        if pairs && saveFiles
            saveas(gcf, strcat(saveFolder, "min vel vs duration correlation.png"))          
        end
    end

end


