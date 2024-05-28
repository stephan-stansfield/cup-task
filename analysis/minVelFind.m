% minVelFind
% 
% Takes a kinematic trajectory and returns the time and value of the first
% velocity local minimum (the inter-peak "dip").

function [timeDip, velDip] = minVelFind(vel,acc,st)
    
    % Velocity local minima and maxima
    velChange = find(diff(sign(acc)))+1;
    
    % Combine pairs of changes where acceleration = 0
    flagInd = [];
    for ind = 1:length(velChange)
        if (acc(velChange(ind)) == 0) && (acc(velChange(ind+1)-1) == 0)
            velChange(ind) = idivide(velChange(ind)+velChange(ind+1),int16(2));
            flagInd = [flagInd,ind];
        end
    end
    velChange(flagInd+1) = [];

    % Count number of velocity changes
    numPeaks = length(velChange);

    % Delete transient changes in velocity from list of peaks & valleys
    velChange = peakFilter(numPeaks, velChange, vel);

    % Recount the number of peaks with short abberations removed
    numPeaks = length(velChange);
    
    % Check that there are at least two velocity peaks
    if numPeaks > 2
        % Return second velocity direction change as minimum velocity, as
        % well as time of occurence
        velDip = vel(velChange(2));
        timeDip = velChange(2)/st;
        
    else
        % If there isn't a local velocity minimum, flag trajectory
        velDip = NaN;
        timeDip = NaN;
    end
    
end