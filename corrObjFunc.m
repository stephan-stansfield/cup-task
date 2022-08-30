% corrObjFunc.m

function [val, gradient] = corrObjFunc(x,plotCorrelation)
   
%     plotCorrelation = true;
    
    B = x(1);
    K = x(2);
    duration_min = 1.00;
    duration_max = 2.00;
    duration_step = 0.01;
    distance = 0.25;
    vStart = 0.15;
    st = 0.001;
    
    correlationArray = nan(length(duration_min:duration_step:duration_max),2);
    i = 1;
    
    % Create model of cart-and-pendulum system coupled with selected impedance
    [sys,Td,zeta] = corrSysCreate(B,K);
    
    for duration = duration_min:duration_step:duration_max
        
        % Simulate motion using generated system and current duration
        [output, dur_corr] = corrSimInputShape(sys,Td,zeta,duration,distance,vStart,st);
        
        vel = output(:,3);
        acc = output(:,5);
        
        % Find inter-peak minimum velocity of simulated profile
        [~, valDip] = minVelFind(vel,acc,st);
        
        
%         % DEBUG
%         disp(strcat('Trial #: ', num2str(i)))

        % Add velocity value and trial duration to array
        correlationArray(i,:) = [dur_corr,valDip];
        
        i = i + 1;
        
    end
    
    % Delete invalid trials from array
    [flagRow, ~] = find(isnan(correlationArray));
    correlationArray(flagRow,:) = [];
    durations = correlationArray(:,1);
    minVels = correlationArray(:,2);

    % Fit linear trendline to data points
%     [poly,S] = polyfit(durations, minVels, 1);
    
    % Calculate Pearson's correlation coefficient, r, and significance
    % level, p
    [R,P] = corrcoef(durations,minVels);
    if length(durations) < 10
    %     if size(R,2) == 1
        disp('Stop here')
        val = -100;
        return
    end
    r = R(1,2);
    p = P(1,2);
    
    if plotCorrelation
        
        labelFontSize = 22;
        axisFontSize = 20;
        titleFontSize = 24;
        
        figure();
        plot(durations,minVels,'o');
        set(gca,'FontSize',axisFontSize)
        xlabel('Trial Duration (s)','FontSize',labelFontSize)
        ylabel('Minimum Velocity (m/s)','FontSize',labelFontSize)
        title({'Interpeak Velocity vs. Duration Correlation', 'Constant B & K'},'FontSize',titleFontSize)
        
    end
    
    % Objective value to be maximized by optimization is correlation
    % coefficient divided by significance. Emphasizes negative correlations
    % that are also significant.
%     val = -r/p;
    if ~isnan(r)
        val = -r;
    else
        val = -100;
    end
    
end