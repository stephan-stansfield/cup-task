% figVAFbars
%
% Create a figure displaying the mean VAF for different simulation types as
% well as error bars

C = readcell('data analysis/VAF summary 2023-03-26.xlsx'); % simulations for 2023 journal paper

% Read all columns:
x = C(1,2:end);                     % x-axis labels
u = 100*cell2mat(C(2,2:end));       % mean VAFs by control model
s = 100*cell2mat(C(3,2:end));       % standard deviations across subjects

X = categorical(x,x);

% Plot settings:
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
% legendFontSize = 20;
marSz = 150;                % Size of scatter plot markers
% markerList = ['o', '*', 'x', 's'];
% % colorList = ['r','b','g'];
% colorList = {[0 0 0], [215 20 20]/256,...
%     [68 114 196]/256, [0.560181 0.691569 0.194885], [41 196 44]/256};
% plotObjects = {};
% 
yellow = [0.9290 0.6940 0.1250];        % MATLAB default yellow
blue = [0, 0.4470, 0.7410];             % MATLAB default blue
green = [0.4660 0.6740 0.1880];         % MATLAB default green

% Assign specific colors and markers to different simulation categories
colorList = {[0 0 0], [0 0 0], [0 0 0], [0 0 0], [0 0 0], [0 0 0]};

markerList = [  "pentagram",... % Original input shaping
                "diamond",...   % Multi-mode input shaping
                "v",...         % Slow mode internal model
                "^",...         % Fast mode internal model
                "square",...    % No impedance internal model (pendulum mode)
                "o"];           % Rigid body internal model (hand impedance mode)

figure();
for i=1:length(X)
    barObjects{i} = errorbar(X(i), u(i), s(i), 'LineStyle','none','LineWidth',2);
    hold on;
    markerObjects{i} = scatter(X(i),u(i),marSz,'Marker',markerList(i));
    % this lines uses the default MATLAB color order for markers:
%     set(markerObjects{i}, 'MarkerFaceColor', get(markerObjects{i},'MarkerEdgeColor')); 
    % this line uses assigned colors from colorList for markers:
    set(markerObjects{i},'MarkerEdgeColor',colorList{i}, 'MarkerFaceColor',colorList{i});
    
    % use default MATLAB color order for bars:
%     set(barObjects{i}, 'Color', get(markerObjects{i},'CData'));
    % use assigned color from colorList:
    set(barObjects{i}, 'Color', colorList{i});
    % put black border on markers:
    set(markerObjects{i}, 'MarkerEdgeColor', 'black', 'LineWidth',1);
end

% % Without explicitly choosing markers & colors:
% for i=1:length(X)
%     plotObjects{i} = scatter(X(i),u(i),100,'filled');
%     % 'LineWidth',2
%     hold on;
%     errorbar(X(i), u(i), s(i),'LineStyle','none','LineWidth',2);
% end

set(gca,'FontSize',axisFontSize)
xlabel('Control Model','FontSize',labelFontSize);
ylabel('Mean Cup Velocity VAF (%)','FontSize',labelFontSize);
ylim([0, 100])
% set(gcf, 'Position',  [100, 100, 700, 600])

% Title
% th = title('Variance Accounted For by Control Model','FontSize',titleFontSize,'FontWeight','Normal');
% titlePos = get(th, 'position');
% titlePos(2) = 107.5;   % Move title up slightly to accomodate adding significance asterisks to figure
% set(th, 'position', titlePos)

% Don't show data labels on x-axis; put in legend instead
set(gca,'box','off')
set(gca, 'xticklabel', []);
set(gca, 'XTick', []);
% % 3 items:
% legend([plotObjects{1} plotObjects{2} plotObjects{3}], x(1:3),...
%     'Location','southeast','FontSize',legendFontSize)

% % 4 items:
% legend([plotObjects{1} plotObjects{2} plotObjects{3} plotObjects{4}], x(1:4),...
%     'Location','southeast','FontSize',legendFontSize)

% 6 items:
% legend([markerObjects{1} markerObjects{2} markerObjects{3} markerObjects{4},...
%     markerObjects{5} markerObjects{6}], x(1:6),...
%     'Location','southeast','FontSize',legendFontSize)
legend([markerObjects{1} markerObjects{2} markerObjects{3} markerObjects{4},...
    markerObjects{5} markerObjects{6}], x(1:6),...
    'Location','southeast')

hold off;


