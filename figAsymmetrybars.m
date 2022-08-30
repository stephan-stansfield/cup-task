% figAsymmetrybars
%
% Create a figure displaying the mean of median peak ratios for different 
% simulation types as well as error bars

C = readcell('best fit simulation/overall peak ratio results.xlsx');

x = C(1,2:end);                 % x-axis labels
u = cell2mat(C(2,2:end));       % mean of median peak ratios
s = cell2mat(C(3,2:end));       % standard deviations of median peak ratios

X = categorical(x,x);

expMean = 1.13;
expStd = 0.16;

% li = expMean*ones(1,14);
% ys = [0:1:13];

% Plot values:
titleFontSize = 28;
labelFontSize = 20;
axisFontSize = 16;
figure();

scatter(X,u);
hold on;
errorbar(X, u, s, 'LineStyle','none');
rectangle('Position',[0,expMean-expStd,15,2*expStd],'FaceColor', [220 220 220]/255,...
    'EdgeColor', [1 1 1])
xlim = get(gca,'xlim');
plot(xlim,[expMean, expMean],'LineStyle','-','Color','black','LineWidth',2);
plot(xlim,[1, 1],'LineStyle','--','Color','black','LineWidth',1.5);
errorbar(X, u, s, 'LineStyle','none','LineWidth',2);
scatter(X,u,'filled');
set(gca,'FontSize',axisFontSize)
title('Blocks 3 & 4','FontSize',titleFontSize);
xlabel('Simulation  Type','FontSize',labelFontSize);
ylabel('Mean of Median Peak Ratios','FontSize',labelFontSize);


