function val = featuresRMSE(acc_RMSE, acc_sim_RMSE, vel_RMSE, vel_sim_RMSE, ...
        theta_RMSE, theta_sim_RMSE, blockStr, blockNum, subjNum, num)
    
    % Load experimental velocity local minima and maxima
%     load(strcat("peak times/", blockStr),'peakArray');
    block3or4 = blockNumFunc(blockNum);
    load(strcat('peak times/',subjNum,'_B',block3or4));
    
    % Get experimental peak timings and values
%     velChange_exp = peakArray{num+1,1}(1,:);
%     velPeaks_exp = peakArray{num+1,1}(2,:);
    velChange_exp = velTimeArray(num+1,:);
    velPeaks_exp = velAmpArray(num+1,:);
    
    % Calculate simulated peak timings and values
    [velChange_sim, ~] = peakCalc(vel_sim_RMSE, acc_sim_RMSE);
    
    % Count the number of peaks with transient changes removed
    numPeaks_exp = length(velChange_exp);
    numPeaks_sim = length(velChange_sim);

%     % Just take the first 3 maxima & minima (two peaks, one valley)
%     if numPeaks_exp > 3
%         velChange_exp = velChange_exp(1:3);
%         numPeaks_exp = 3;
%     end
%     if numPeaks_sim > 3
%         velChange_sim = velChange_sim(1:3);
%         numPeaks_sim = 3;
%     end
    
    if numPeaks_sim > numPeaks_exp
        velChange_sim = velChange_sim(1:numPeaks_exp);
    end

%     % Alternative: *only* take the inter-peak velocity minimum
%     if numPeaks_exp >= 2
%         velChange_exp = velChange_exp(2);
%         numPeaks_exp = 1;
%     end
%     if numPeaks_sim >= 2
%         velChange_sim = velChange_sim(2);
%         numPeaks_sim = 1;
%     end
        
    % Reject profiles that don't have the same number of peaks as the 
    % experimental trial
    if numPeaks_exp ~= numPeaks_sim
        val = 100000;

%         velPeaks_exp = vel_RMSE(velChange_exp);                             % comment
        velPeaks_sim = vel_sim_RMSE(velChange_sim);

    else

        % Velocity values at the local minima & maxima
%         velPeaks_exp = vel_RMSE(velChange_exp);                             % comment
        velPeaks_sim = vel_sim_RMSE(velChange_sim);

        % Residual angle
        res_exp         = theta_RMSE(end);
        res_sim         = theta_sim_RMSE(end);

        % Trial duration
        dur_exp         = length(vel_RMSE);
        dur_sim         = length(vel_sim_RMSE);

        % Calculate difference between experimental and simulated values
        diffPeaks       = velPeaks_exp' - velPeaks_sim;
        diffPeakTimes   = velChange_exp' - velChange_sim*0.001;             % DEBUG: 0.001 is st here. Have to change if changing simulation step size!
        diffRes         = res_exp - res_sim;
        diffDur         = dur_exp - dur_sim;   

        % Normalize differences to range of experimental trajectory
        diffPeaks       = diffPeaks/(max(vel_RMSE)-min(vel_RMSE));
        diffPeakTimes   = diffPeakTimes/length(vel_RMSE);
        diffRes         = diffRes/(max(theta_RMSE)-min(theta_RMSE));
        diffDur         = diffDur/max(dur_exp,dur_sim);

        % Square the difference
        sqPeaks         = diffPeaks.^2;
        sqPeakTimes     = diffPeakTimes.^2;
        sqRes           = diffRes.^2;
        sqDur           = diffDur.^2;

        % Normalize velocity values and times by number of peaks
        % (so that residual angle is just as important as the
        % velocity values and times)
        sqPeaks         = sqPeaks/numPeaks_exp;
        sqPeakTimes     = sqPeakTimes/numPeaks_exp;

        % Assign weights to each feature
        wPeaks          = 1.0;
        wPeakTimes      = 1.0;
        wRes            = 0;
        wDur            = 0;

        % Compile terms into one vector and apply variable weights
        diffVector      = vertcat(wPeaks*sqPeaks,wPeakTimes*sqPeakTimes,...
            wRes*sqRes,wDur*sqDur);

        % Take square root of sum of mean squared errors
        val = sqrt(1/length(diffVector)*sum(diffVector));
                
    end

    % Store location of peaks in globalData file
    globalData('velChange_exp', velChange_exp);
    globalData('velPeaks_exp', velPeaks_exp);
    globalData('velChange_sim', velChange_sim);
    globalData('velPeaks_sim', velPeaks_sim);
        
end