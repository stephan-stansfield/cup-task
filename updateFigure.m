function updateFigure(figHand, d)
% figHand.ZData(d(1)) = d(2);
% drawnow('limitrate');

set(groot,'CurrentFigure',figHand);
subplot(plotRows,plotCols,num-numStart+1)                               
plot(te,pos,'LineWidth',2)
hold on;
if plotSim
    plot(tc,pos_sim,'LineWidth',2)
end
trial = num2str(num+1);
title(['Trial' ' ' trial]);
xlabel('t (s)')
ylabel('cart position (m)')

end