% Calculates the timings and amplitudes of velocity local minima and
% maxima. Calls peakFilter to remove transient changes in velocity (changes
% that happen under some minimum time threshold).
function [velChangeTimesSim, velChangeAmpsSim] = peakCalc(vel_sim, acc_sim)

    % Find velocity local minima and maxima
    velChangeTimesSim = find(diff(sign(acc_sim)))+1;
%     velChange_exp = find(diff(sign(acc_RMSE)))+1;

    % Count number of velocity changes
    numPeaks_sim = length(velChangeTimesSim);
%     numPeaks_exp = length(velChange_exp);

    % Delete transient changes in velocity from list of peaks & valleys
    velChangeTimesSim = peakFilter(numPeaks_sim, velChangeTimesSim, vel_sim);
%     velChange_exp = peakFilter(numPeaks_exp, velChange_exp, vel_RMSE);
      
    % Velocity values at the local minima & maxima
    velChangeAmpsSim = vel_sim(velChangeTimesSim);
%     velPeaks_exp = vel_RMSE(velChange_exp);
        
end