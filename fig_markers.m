% figVAFbars
%
% Create a figure displaying the mean VAF for different simulation types as
% well as error bars

C = readcell('data analysis/VAF summary 2022-10-17.xlsx');                  % simulations for SfN 2022 poster

% Read all columns:
x = C(1,2:end);                     % x-axis labels
u = 100*cell2mat(C(2,2:end));       % means of median VAFs
s = 100*cell2mat(C(3,2:end));       % standard deviations

% Show just original, multi-mode, and simplified internal model:
% x = C(1,2:4);                     % x-axis labels
% u = 100*cell2mat(C(2,2:4));       % means (of median VAFs)
% s = 100*cell2mat(C(3,2:4));       % standard deviations

% Show the 3 input shaping models + submovement fits:
% x = C(1,2:5);                     % x-axis labels
% u = 100*cell2mat(C(2,2:5));       % means (of median VAFs)
% s = 100*cell2mat(C(3,2:5));       % standard deviations

X = categorical(x,x);

% Plot settings:
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;
marSz = 500;                % Size of scatter plot markers
% markerList = ['o', '*', 'x', 's'];
% % colorList = ['r','b','g'];
% colorList = {[0 0 0], [215 20 20]/256,...
%     [68 114 196]/256, [0.560181 0.691569 0.194885], [41 196 44]/256};
% plotObjects = {};
% 
orange = [0.9290 0.6940 0.1250];        % MATLAB default orange
green = [0.4660 0.6740 0.1880];         % MATLAB default green

% Assign specific colors and markers to different simulation categories
colorList = {   [0 0 0],...                                     % Original input shaping
                orange, orange, orange, orange, orange,...     % No feedforward
                green, green, green, green, green};     % Feedforward

markerList = [  "x",...                         % Human subject data
                "pentagram",...                 % Original input shaping
                "diamond",...                   % Multi-mode input shaping
                'v',...                         % Slow mode internal model
                '^',...                         % Fast mode internal model
                'o',...                         % No impedance internal model (pendulum mode)
                "square"];                      % Rigid body internal model (hand impedance mode)
%                 'v', '^', 'o', "square"];       % Repeat for feedforward

figure();
% First plot markers and legend
for i=1:length(markerList)
    plotObjects{i} = scatter(X(i),u(i),marSz,'Marker',markerList(i),...
        'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0],'LineWidth',2)
    hold on;
end

% % Without explicitly choosing markers & colors:
% for i=1:length(X)
%     plotObjects{i} = scatter(X(i),u(i),100,'filled');
%     % 'LineWidth',2
%     hold on;
%     errorbar(X(i), u(i), s(i),'LineStyle','none','LineWidth',2);
% end

set(gca,'FontSize',axisFontSize)
xlabel('Simulation  Type','FontSize',labelFontSize);
ylabel('Mean of Median VAF (%)','FontSize',labelFontSize);
% ylim([0, 102.5])
% set(gcf, 'Position',  [100, 100, 700, 600])

th = title('Variance Accounted For by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal');
titlePos = get(th, 'position');
titlePos(2) = 105;   % Move title up slightly to accomodate adding significance askterisks to figure
set(th, 'position', titlePos)

% % Don't show data labels on x-axis; put in legend instead
set(gca,'xticklabel',[])
% % 3 items:
% legend([plotObjects{1} plotObjects{2} plotObjects{3}], x(1:3),...
%     'Location','southeast','FontSize',legendFontSize)

% 4 items:
% legend([plotObjects{1} plotObjects{2} plotObjects{3} plotObjects{4}], x(1:4),...
%     'Location','southeast','FontSize',legendFontSize)

hold off;

