function velChangeTimes = peakFilter(numPeaks, velChangeTimes, vel)
    
    % Delete transient changes in velocity, defined as less than
    % 100 ms between adjacent acceleration sign changes.
    for indPeak = 1:(numPeaks-2)

        % If there are 3 adjacent direction changes with < 100 ms between
        % them, flag the middle one for deletion. Then, flag the one that
        % has a smaller space between it and the next one.
        if velChangeTimes(indPeak+1)-velChangeTimes(indPeak) < 100 &&...
                velChangeTimes(indPeak+2)-velChangeTimes(indPeak+1) < 100

            % Flag the middle time change from the array
            velChangeTimes(indPeak+1) = -100;

            % If the velocity change is a local maximum, flag the smaller
            % value for deletion. If a local minimum, flag the larger value
            % for deletion. Index of middle direction change will be even
            % for local maxima and odd for local minima. (This is because
            % direction changes always alternate and velocity always starts
            % in positive direction.)
            if mod(indPeak+1,2) == 0 % Local maximum
                % Flag smaller value
                if vel(velChangeTimes(indPeak)) < vel(velChangeTimes(indPeak+2))
                    velChangeTimes(indPeak) = -100;
                else
                    velChangeTimes(indPeak+2) = -100;
                end
            else % Local minimum
                % Flag larger value
                if vel(velChangeTimes(indPeak)) < vel(velChangeTimes(indPeak+2))
                    velChangeTimes(indPeak+2) = -100;
                else
                    velChangeTimes(indPeak) = -100;
                end
            end
        end
    end

    % Delete flagged indices from array and recalculate number of peaks
    velChangeTimes = velChangeTimes(velChangeTimes ~= -100);
    numPeaks = length(velChangeTimes);
    
    % If there are only 2 direction changes separated by < 100 ms, flag
    % both for deletion
    for indPeak = 1:(numPeaks-1)
        if velChangeTimes(indPeak+1)-velChangeTimes(indPeak) < 100
            velChangeTimes(indPeak) = -100;
            velChangeTimes(indPeak+1) = -100;
        end
    end
    
    % Delete flagged indices from array
    velChangeTimes = velChangeTimes(velChangeTimes ~= -100);