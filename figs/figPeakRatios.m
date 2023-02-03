% figPeakRatios
%
% Create a figure displaying the mean of median peak ratios for different 
% simulation types as well as error bars

% C = readcell('data analysis/peak ratio results overall.xlsx');
% C = readcell('data analysis/peak ratio results FF only 2022-11-03.xlsx');
C = readcell('data analysis/peak ratio results select 2022-11-22.xlsx');

x = C(1,2:end);                 % x-axis labels
u = cell2mat(C(2,2:end));       % mean of median peak ratios
s = cell2mat(C(3,2:end));       % standard deviations of median peak ratios

X = categorical(x,x);

expMean = 1.13;
expStd = 0.16;

% li = expMean*ones(1,14);
% ys = [0:1:13];

% Plot settings
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;
marSz = 250;                    % Size of scatter plot markers
green = [0.4660 0.6740 0.1880]; % MATLAB default green
blue = [0, 0.4470, 0.7410];     % MATLAB default blue
markerList = [  "x",...         % Human subject data
                "pentagram",... % Uncoupled system (original input shaping)
                "v",...         % Slow mode, NF
                "o",...         % No impedance, NF
                "diamond",...   % Multi-mode, FF
                "v",...         % Slow mode, FF
                "square",...    % Rigid body, FF
                "o"];           % No impedance, FF

colorList = {   [0 0 0], [0 0 0],...            % Human data, uncoupled system
                blue, blue,...                  % No feedforward
                green, green, green, green};    % Feedforward

figure();
% barObjects{1} = errorbar(X(1), u(1), s(1),'LineStyle','none','LineWidth',2);
% hold on;
% markerObjects{1} = scatter(X(1),u(1),200,'Marker',markerList(1),...
%     'LineWidth',1.5);
% set(markerObjects{1}, 'MarkerFaceColor', get(markerObjects{1},'MarkerEdgeColor'));
% % set(barObjects{1}, 'Color', get(markerObjects{1},'CData'));
% set(barObjects{1}, 'Color', 'black');
% set(markerObjects{1}, 'MarkerEdgeColor', 'black', 'LineWidth',1);

%         'MarkerEdgeColor',"k",'MarkerFaceColor',"k",'LineWidth',1.5);


for i = 1:length(X)
    barObjects{i} = errorbar(X(i), u(i), s(i), 'LineStyle','none','LineWidth',2);
    hold on;
    markerObjects{i} = scatter(X(i),u(i),marSz,'Marker',markerList(i))
        set(markerObjects{i}, 'MarkerFaceColor', colorList{i},'MarkerEdgeColor', colorList{i})
%     set(markerObjects{i}, 'MarkerFaceColor', get(markerObjects{i},'MarkerEdgeColor'));
%         'MarkerEdgeColor',green,'MarkerFaceColor',green);
%     hold on;
    
%     set(barObjects{i}, 'Color', get(markerObjects{i},'CData'));
    set(barObjects{i}, 'Color', colorList{i});
    set(markerObjects{i}, 'MarkerEdgeColor', 'black', 'LineWidth',1);
%     errorbar(X(i), u(i), s(i), 'Color',green,'LineStyle','none','LineWidth',2);
end
% hold on;
% errorbar(X, u, s, 'Color',green,'LineStyle','none');

rectX = [0.75, 8.25]
plot(rectX,[1, 1],'LineStyle','--','Color','black','LineWidth',1.5);          % 1 (no asymmetry) as dashed line
plot(rectX,[expMean, expMean],'LineStyle','-','Color','black','LineWidth',1.5); % experimental mean of medians as solid line
a = gca
ch = get(gca,'Children')
set(gca,'Children',[ch(17) ch(18) ch(15) ch(16) ch(13) ch(14) ch(11) ch(12) ch(9) ch(10) ch(7) ch(8)...
    ch(5) ch(6) ch(3) ch(4) ch(1) ch(2) ]);
ch2 = get(gca,'Children')

% axPos = get(gca,'Position')
% xlim = get(gca,'xlim');

% rectangle('Position',[0,expMean-expStd,15,2*expStd],'FaceColor', [220 220 220]/255,...
%     'EdgeColor', [1 1 1])
% rectangle('Position',[xlim(1),expMean-expStd,axPos(3),2*expStd],'FaceColor', [220 220 220]/255,...
%     'EdgeColor', [1 1 1])
% rectangle('Position',[0,expMean-expStd,5,2*expStd],'FaceColor', [220 220 220]/255,...
%     'EdgeColor', [1 1 1])
% rectX = xlim

% rectX = [0.75, 6.25]

% rectY = expMean + expStd*[1 -1]
% patch(rectX([1,2,2,1]), rectY([1,1,2,2]),'k', ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.2);

% plot(rectX,[1, 1],'LineStyle','--','Color','black','LineWidth',1.5);          % 1 (no asymmetry) as dashed line

% plot(rectX,[expMean, expMean],'LineStyle','-','Color','black','LineWidth',2); % experimental mean of medians as solid line
% errorbar(X, u, s,'Color',green,'LineStyle','none','LineWidth',2);
% hold on;
% for i = 1:length(X)
%     scatter(X(i),u(i),marSz,'Marker',markerList(i),...
%         'MarkerEdgeColor',green,'MarkerFaceColor',green);
% end
% scatter(X,u,'filled');
set(gca,'FontSize',axisFontSize)
title('Peak Ratio by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal');
xlabel('Simulation  Type','FontSize',labelFontSize);
ylabel('Mean of Median Peak Ratios','FontSize',labelFontSize);
set(gca,'xticklabel',[]) % Don't show data labels on x-axis; put in legend instead

% Add legend and set desired marker size
% legend([markerObjects{1} markerObjects{2} markerObjects{3} markerObjects{4},...
%     markerObjects{5} markerObjects{6}], x(1:6),...
%     'Location','southoutside','FontSize',legendFontSize)
% [~, objh] = legend([markerObjects{1} markerObjects{2} markerObjects{3} markerObjects{4},...
%     markerObjects{5} markerObjects{6}], x(1:6),...
%     'Location','southoutside','FontSize',legendFontSize)
% objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
% set(objhl, 'Markersize', 100); %// set marker size as desired



