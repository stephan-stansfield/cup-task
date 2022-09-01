function updateFigure(figHand, d)
% figHand.ZData(d(1)) = d(2);
% drawnow('limitrate');

% DEBUG
% disp('Figure handle: ')
% disp(figHand)

% Unpack sent data
plotRows = d(1);
plotCols = d(2);
num = d(3);
numStart = d(4);
len_te = d(5);
len_tc = d(6);
te = d(7:len_te+6);
pos = d(len_te+7:2*len_te+6);
tc = d(2*len_te+7:2*len_te+len_tc+6);
pos_sim = d(2*len_te+len_tc+7:end);

% DEBUG
% disp('Size of te:')
% size(te)
% disp('Size of pos:')
% size(pos)
% 
% disp('Size of tc:')
% size(tc)
% disp('Size of pos_sim:')
% size(pos_sim)

% Set current figure to the specified handle and plot
set(groot,'CurrentFigure',figHand);
subplot(plotRows,plotCols,num-numStart+1)                               
plot(te,pos,'LineWidth',2)
hold on;
plot(tc,pos_sim,'LineWidth',2)

trial = num2str(num+1);
title(['Trial' ' ' trial]);
xlabel('t (s)')
ylabel('cart position (m)')

end