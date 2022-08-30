% figPlotAllTrials
%
% Produces a figure showing individual trials in one block as subplots.
% This is the "standard" figure that I have been using to examine the
% results of the optimization algorithm.
%
% This code was originally built into optimization.m. I put it in a
% separate script so that it can be called by other functions, e.g., to
% plot previously simulated best-fit trials including velocity peak markers
% found by peakAsymmetry.m

function figPlotAllTrials(cv, plotRows, plotCols, num, numStart, te, vel,...
        plotSim, plotExp, plotDes, tc, velSim, plotPeaks, velChangeTimesSim,...
        velChangeTimesExp, st, velChangeAmpsSim, velChangeAmpsExp, velDes)
    % Cart velocity (vel): plot experimental versus simulation
        set(groot,'CurrentFigure',cv);
        subplot(plotRows,plotCols,num-numStart+1)   
        if plotExp
            plot(te,vel,'LineWidth',2)
            hold on;
            
            if plotPeaks && sum(isnan(velChangeTimesExp)) == 0
%                 plot(velChangeTimesExp*st,velChangeAmpsExp,'Marker','s',...
%                     'MarkerSize',12,'Color','magenta','LineWidth',1.5,...
%                     'LineStyle','none')
%                 plot(te(velChangeTimesExp),velChangeAmpsExp,'Marker','s',...
%                     'MarkerSize',12,'Color','magenta','LineWidth',1.5,...
%                     'LineStyle','none')
                plot(velChangeTimesExp,velChangeAmpsExp,'Marker','s',...
                    'MarkerSize',12,'Color','magenta','LineWidth',1.5,...
                    'LineStyle','none')
            end
        end
        
        if plotSim && ~plotDes
            plot(tc,velSim,'LineWidth',2)
            hold on;

            % Plot simulated local velocity max & min
            if plotPeaks && sum(isnan(velChangeTimesSim)) == 0  
                plot(velChangeTimesSim*st,velChangeAmpsSim,'Marker','s',...
                    'MarkerSize',12,'Color','blue','LineWidth',1.5,...
                    'LineStyle','none')
            end
        end

        if plotDes
            % Plot individual input velocity submovements
            [vcRows, vcCols] = size(velDes);
            td = 0:st:(vcCols-1)*st;
            for row = 1:vcRows
                plot(td,velDes(row,:),'LineWidth',1.25,'Color',[0.9290, 0.6940, 0.1250])   
            end
            % Plot simulated output velocity
            plot(tc, velSim, 'LineWidth', 2,'LineStyle','--','Color',[0.8500, 0.3250, 0.0980]);
        end
        trial = num2str(num);
        title(['Trial' ' ' trial]);
        xlabel('Time (s)')
        ylabel('Cart Velocity (m/s)')