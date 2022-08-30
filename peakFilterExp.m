function velChangeTimes = peakFilterExp(numPeaks, velChangeTimes, vel, t)
    
%     % Initialize array to hold indices of markers to be deleted
%     flagTimes = [];

    % Initialize array of even & odd numbers to indicate peaks & valleys
    peakOrValley = zeros(1,numPeaks);
    for i = 1:numPeaks
       if mod(i,2) == 1 % odd
           peakOrValley(i) = 1;
       else
           peakOrValley(i) = -1;
       end
    end

    % Delete transient changes in velocity, defined as less than
    % 100 ms between adjacent acceleration sign changes.
    for indPeak = 1:(numPeaks-2)

        % DEBUG
%         disp('Peak index: ');
%         disp(indPeak)
        
        % If there are 3 adjacent direction changes with < 100 ms
        % between them, flag the middle one for deletion. Then,
        % flag the one that has a smaller space between it and the 
        % next one.
        if t(velChangeTimes(indPeak+1)) - t(velChangeTimes(indPeak)) < 0.100 &&...
                t(velChangeTimes(indPeak+2)) - t(velChangeTimes(indPeak+1)) < 0.100
            
            % DEBUG
%             disp('3 adjacent peaks!')

            % Flag the middle time change from the array
%             velChangeTimes(foo+1) = -100;
%             flagTimes = [flagTimes, velChangeTimes(indPeak+1)];
            % Delete the middle time change from the array
            velChangeTimes(indPeak+1) = [];
            peakOrValley(indPeak+1) = [];

            % If the velocity change is a local maximum, flag the
            % smaller value for deletion. If a local minimum, flag
            % the larger value for deletion.
            % Index of middle direction change will be even for
            % local maxima and odd for local minima. (This is
            % because direction changes always alternate and
            % velocity always starts in positive direction.)
            if peakOrValley(indPeak+1) == 1    % Local maximum
                % Flag smaller value
                if vel(velChangeTimes(indPeak)) < vel(velChangeTimes(indPeak+1))
%                     velChangeTimes(foo) = -100;
%                     flagTimes = [flagTimes, velChangeTimes(indPeak)];
                    velChangeTimes(indPeak) = [];
                    peakOrValley(indPeak) = [];
                else
%                     velChangeTimes(foo+2) = -100;
%                     flagTimes = [flagTimes, velChangeTimes(indPeak+2)];
                    velChangeTimes(indPeak+1) = []; % Note that this was indPeak+2, but we deleted the value at indPeak+1
                    peakOrValley(indPeak+1) = [];
                end
            else    % Local minimum
                % Flag larger value
                if vel(velChangeTimes(indPeak)) < vel(velChangeTimes(indPeak+1))
%                     velChangeTimes(foo+2) = -100;
%                     flagTimes = [flagTimes, velChangeTimes(indPeak+2)];
                    velChangeTimes(indPeak+1) = []; % Note that this was indPeak+2, but we deleted the value at indPeak+1
                    peakOrValley(indPeak+1) = [];
                else
%                     velChangeTimes(foo) = -100;
%                     flagTimes = [flagTimes, velChangeTimes(indPeak)];
                    velChangeTimes(indPeak) = [];
                    peakOrValley(indPeak) = [];
                end
            end
        end
        
        % Check if loop has reached end of re-sized array
%         % DEBUG
%         disp('Peak index after trimming: ')
%         disp(indPeak)
%         disp('Length of velChangeTimes: ')
        if indPeak >= length(velChangeTimes)-2
            break
        end
        
    end

    % Delete flagged indices from array and recalculate number of peaks
%     velChangeTimes = velChangeTimes(velChangeTimes ~= -100);
%     velChangeTimes = setdiff(velChangeTimes,flagTimes);
    
    numPeaks = length(velChangeTimes);
    
    % If there are only 2 direction changes separated by < 100
    % ms, flag them both for deletion
    for indPeak = 1:(numPeaks-1)
        if t(velChangeTimes(indPeak+1)) - t(velChangeTimes(indPeak)) < 0.100
%             velChangeTimes(indPeak) = -100;
%             velChangeTimes(indPeak+1) = -100;
%             flagTimes = [flagTimes, velChangeTimes(indPeak)];
%             flagTimes = [flagTimes, velChangeTimes(indPeak+1)];
            velChangeTimes(indPeak+1) = [];
            velChangeTimes(indPeak) = [];
            peakOrValley(indPeak+1) = [];
            peakOrValley(indPeak) = [];
        end
        
        % Check if loop has reached end of re-sized array
        if indPeak >= length(velChangeTimes)-1
            break
        end
        
    end
    
    % Delete flagged indices from array
%     velChangeTimes = velChangeTimes(velChangeTimes ~= -100);
%     velChangeTimes = setdiff(velChangeTimes,flagTimes);
