% figVAFbars
%
% Create a figure displaying the mean VAF for different simulation types as
% well as error bars

% C = readcell('best fit simulation/overall VAFs.xlsx');                        % file containing all 11 simulation types
% C = readcell('best fit simulation/overall VAFs select RQE.xlsx');             % file containing only 4 simulation types
% C = readcell('best fit simulation/overall VAFs select RQE std of means.xlsx'); % file containing only 4 simulation types taking standard deviation of subject mean VAFs
% C = readcell('data analysis/VAF summary all models.xlsx');                  % include all models
C = readcell('data analysis/VAF summary 2022-08-17.xlsx');                  % only show 3 or 4 simulations


% Read all columns:
% x = C(1,2:end);                     % x-axis labels
% u = 100*cell2mat(C(2,2:end));       % means (of median VAFs)
% s = 100*cell2mat(C(3,2:end));       % standard deviations

% Show just original, multi-mode, and simplified internal model:
x = C(1,2:4);                     % x-axis labels
u = 100*cell2mat(C(2,2:4));       % means (of median VAFs)
s = 100*cell2mat(C(3,2:4));       % standard deviations

% Show the 3 input shaping models + submovement fits:
% x = C(1,2:5);                     % x-axis labels
% u = 100*cell2mat(C(2,2:5));       % means (of median VAFs)
% s = 100*cell2mat(C(3,2:5));       % standard deviations

X = categorical(x,x);

% Plot values:
titleFontSize = 26;
labelFontSize = 24;
axisFontSize = 16;
legendFontSize = 20;
markerList = ['o', '*', 'x', 's'];
% colorList = ['r','b','g'];
colorList = {[0 0 0], [215 20 20]/256,...
    [68 114 196]/256, [0.560181 0.691569 0.194885], [41 196 44]/256};
plotObjects = {};

figure();
% First plot markers and legend
for i=1:length(X)
    plotObjects{i} = scatter(X(i),u(i),100,'filled','Marker',markerList(i),...
        'MarkerEdgeColor',colorList{i},'MarkerFaceColor',colorList{i});
    % 'LineWidth',2
    hold on;
    errorbar(X(i), u(i), s(i),'Color',colorList{i},'LineStyle','none','LineWidth',2);
end

% Without explicitly choosing markers & colors:
% for i=1:length(X)
%     plotObjects{i} = scatter(X(i),u(i),100,'filled');
%     % 'LineWidth',2
%     hold on;
%     errorbar(X(i), u(i), s(i),'LineStyle','none','LineWidth',2);
% end

set(gca,'FontSize',axisFontSize)
xlabel('Simulation  Type','FontSize',labelFontSize);
ylabel('Mean of Median VAF (%)','FontSize',labelFontSize);
ylim([0, 102.5])
set(gcf, 'Position',  [100, 100, 700, 600])

th = title('Variance Accounted for by Simulation Type','FontSize',titleFontSize,'FontWeight','Normal');
titlePos = get(th, 'position');
titlePos(2) = 105;   % Move title up slightly to accomodate adding significance askterisks to figure
set(th, 'position', titlePos)

% Don't show data labels on x-axis; put in legend instead
set(gca,'xticklabel',[])
% 3 items:
legend([plotObjects{1} plotObjects{2} plotObjects{3}], x(1:3),...
    'Location','southeast','FontSize',legendFontSize)

% 4 items:
% legend([plotObjects{1} plotObjects{2} plotObjects{3} plotObjects{4}], x(1:4),...
%     'Location','southeast','FontSize',legendFontSize)

hold off;



%%%%
% % Original Input Shaping
% x1 = 'IS';      % x-axis label
% u1 = 0.253;     % mean of all trials
% s1 = 0.581;     % standard deviation of all trials
%
% % Multi-Mode Input Shaping, No Feedforward
% x2 = 'MM-NF';
% u2 = 0.141;
% s2 = 0.789;
% 
% % No-Impedance Simplification, No Feedforward
% x3 = 'NI-NF';
% u3 = 0.187;
% s3 = 0.645;
% 
% % Rigid-Body Simplification, No Feedforward
% x4 = 'RB-NF';
% u4 = 0.022;
% s4 = 0.617;
% 
% % Slow Mode Simplification, No Feedforward
% x5 = 'SM-NF';
% u5 = 0.433;
% s5 = 0.542;
% 
% % Fast Mode Simplification, No Feedforward
% x6 = 'FM-NF';
% u6 = -0.26;
% s6 = 0.93;
% 
% % Multi-Mode Input Shaping, Feedforward
% x7 = 'MM-F';
% u7 = 0.502;
% s7 = 0.439;
% 
% % No-Impedance Simplification, Feedforward
% x8 = 'NI-F';
% u8 = 0.614;
% s8 = 0.371;
% 
% % Rigid-Body Simplification, Feedforward
% x9 = 'RB-F';
% u9 = 0.672;
% s9 = 0.220;
% 
% % Slow Mode Simplification, Feedforward
% x10 = 'SM-F';
% u10 = 0.620;
% s10 = 0.443;
% 
% % Fast Mode Simplification, Feedforward
% x11 = 'FM-F';
% u11 = 0.665;
% s11 = 0.284;
% 
% % Submovement Fitting
% x12 = 'SB';
% u12 = 0.919;
% s12 = 0.078;
% 
% x = {x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12};
% u = 100*[u1; u2; u3; u4; u5; u6; u7; u8; u9; u10; u11; u12];
% s = 100*[s1; s2; s3; s4; s5; s6; s7; s8; s9; s10; s11; s12];


