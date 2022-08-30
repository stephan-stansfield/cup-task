% Trims experimental trials to proper start and end indices and saves to a
% new file along with the average of each kinematic variable.

function statisticsExpMean()
    
numStart    = 1;
numEnd      = 50;
st          = 0.001;
    
for blockNum = 1:22

    % Get information about experimental block
    [~, subjStr, trialDate, trialStr, blockStr, ~,...
        ~, ~] = blockDictionary(blockNum);
    
    % Load trimming indices
    load(blockStr,'Expression1')
    
    saveFolder = strcat("experimental data trimmed/",subjStr,trialDate,trialStr);
    saveFolder = erase(saveFolder, "_trial_");
    
    if(not(isfolder(saveFolder)))
        mkdir(saveFolder);
    end
    
    % Loop through trials, trimming each one and saving new values and
    % averages in folder
    for num = numStart:1:numEnd
       
        numStr = num2str(num-1);

        % Construct file name and load file of one experimental trial
        fileStr = strcat(subjStr,trialDate,trialStr,numStr);
        load(fileStr,'pos','theta','vel','omega','acc','alpha','t')
        
        numStr = num2str(num);
        
        % Construct file name to save to
        saveFile = strcat(saveFolder,"/",numStr);
        
        % Load start and stop indices of corresponding trial
        start   = Expression1(num,1);
        stop    = Expression1(num,2);

        % Trim experimental data
        [pos,theta,vel,omega,acc,alpha,tdes] = trimData(pos,theta,vel,omega,...
            acc,alpha,t,st,start,stop);
        
        % Take average values of each kinematic variable
        posAvg      = mean(pos);
        thetaAvg    = mean(theta);
        velAvg      = mean(vel);
        omegaAvg    = mean(omega);
        accAvg      = mean(acc);
        alphaAvg    = mean(alpha);
        
        % Save trimmed trajectories and average values to a new file
        save(saveFile,'pos','theta','vel','omega','acc','alpha','tdes',...
            'posAvg','thetaAvg','velAvg','omegaAvg','accAvg','alphaAvg')
        
    end
    
end