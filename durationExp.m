% durationExp

blockStart = 1;
blockEnd = 22;
trialStart = 0;
trialEnd = 49;

% Initialize array to hold all trial durations
durationArray = nan(1100,1);

for  blockNum = blockStart:1:blockEnd
    
    % Get info for current block
    [subjNum, subjStr, trialDate, trialStr, blockStr, ~, ~,invalidTrials] = ...
        blockDictionary(blockNum);
    
    % Load start and stop indices of trials for this block
    load(strcat("trim times/",blockStr),'Expression1')

    for trialNum = trialStart:1:trialEnd
        
        % Load file of one experimental trial
        numStr = num2str(trialNum);
        fileStr = strcat(subjStr,trialDate,trialStr,numStr);
        load(fileStr,'pos','theta','vel','omega','acc','alpha','t');
        
        % Load start and stop indices of corresponding trial. Note that
        % trials are indexed from 0, so add 1 to access correct row in array
        start   = Expression1(trialNum+1,1);
        stop    = Expression1(trialNum+1,2);

        % Calculate duration of each trimmed experimental trial
        trialDur = t(stop) - t(start);
        
        % Add duration to array
        durationArray((blockNum-1)*50+(trialNum+1)) = trialDur;
        
    end

end