% statisticsVAF
%
% Calculates Variance Accounted For (VAF) using same equation as
% Svinin et al (2019)

clear all;
close all;
clc;

tic % Time file execution

addpath('data', 'trim times', genpath('experimental data'));
st          = 0.001;                                                        % Simulation time step (s)
stNorm      = 0.0001;                                                       % Normalized time step (s)
tNorm       = 0:stNorm:1;                                                   % Normalized time vector
normLen     = length(tNorm);    

% Choose whether to take average of all experimental and simulated
% trials in block, or to take early and late simulated trials only.
type = "all";
%     type = "early and late";

% Choose whether to only calculate velocity VAF, or VAF for all kinematic
% variables
justVel = true;

% Choose whether to save generated statistics
saveFile = true;

% Choose which experimental blocks to analyze (all acceptable blocks
% would be 1-22)
blockStart  = 1;
blockEnd    = 22;
numBlocks   = blockEnd - blockStart + 1;

% Choose range of experimental trials to analyze
numStart    = 1;
numEnd      = 50;

% Choose range of simulation types to analyze
simStart = 1;
simEnd = 6;

% Label for simulation type
simLabels = {'IS', 'MM-FF', 'SM-FF', 'FM-FF', 'RB-FF', 'NI-FF'};
   

% Create table to hold overall VAF mean and standard deviation values for
% all simulations analyzed
overallResults = {'Simulation Type'; 'Mean VAF'; 'Standard Deviation of Mean VAFs';...
    'Mean of Median VAFs'; 'Std Dev of Median VAFs'};

% Create arrays to hold mean and std dev VAF data for all simulations
meanVAF = nan(simEnd, blockEnd/2);
sdevVAF = nan(simEnd, blockEnd/2);
medianVAF = nan(simEnd, blockEnd/2);
sdevMedianVAF = nan(simEnd, 1);
meanOfMedianVAF = nan(simEnd, 1);

% Loop through simulation types
for simNum = simStart:simEnd
    
    % 1. Nominal Input Shaping (No Impedance Model)
    % 2. Multi-Mode Input Shaping, no FF
    % 3. Slow Mode IS, no FF
    % 4. Fast Mode IS, no FF
    % 5. Rigid Body IS, no FF
    % 6. No Impedance IS, no FF
    % 7. Multi-Mode Input Shaping, FF
    % 8. Slow Mode IS, FF
    % 9. Fast Mode IS, FF
    % 10. Rigid Body IS, FF
    % 11. No Impedance IS, FF
    
    % Get name of folder holding best-fit simulations for selected
    % simulation type
    parentFolder = simulationFolder(simNum);
    if parentFolder == "/"
        continue
    end
    
    % Initialize arrays to hold VAF values
%     VAFVel = zeros(50,numBlocks);
%     blockAvgVAFVel = zeros(1,numBlocks);
%     subjAvgVAFVel = zeros(1,numBlocks);
    VAFVel = nan(50,numBlocks);
    blockAvgVAFVel = nan(1,numBlocks);
    subjAvgVAFVel = nan(1,numBlocks);
    simulationAvgVAFVel = nan(1,numBlocks);
    subjSDVel = nan(1,numBlocks);
    subjAvgSDVel = nan(1,numBlocks);
    simulationSDVel = nan(1,numBlocks);
    subjMedianVAF = nan(1,numBlocks);
    subjMedianSD = nan(1,numBlocks);
    
    if ~justVel
        VAFPos = zeros(50,numBlocks);
        VAFAcc = zeros(50,numBlocks);
        VAFTheta = zeros(50,numBlocks);
        VAFOmega = zeros(50,numBlocks);
        VAFAlpha = zeros(50,numBlocks);

        blockAvgVAFPos = zeros(1,numBlocks);
        blockAvgVAFAcc = zeros(1,numBlocks);
        blockAvgVAFTheta = zeros(1,numBlocks);
        blockAvgVAFOmega = zeros(1,numBlocks);
        blockAvgVAFAlpha = zeros(1,numBlocks);

        subjAvgVAFPos = zeros(1,numBlocks);
        subjAvgVAFAcc = zeros(1,numBlocks);
        subjAvgVAFTheta = zeros(1,numBlocks);
        subjAvgVAFOmega = zeros(1,numBlocks);
        subjAvgVAFAlpha = zeros(1,numBlocks);

        simulationAvgVAFPos = nan(1,numBlocks);
        simulationAvgVAFAcc = nan(1,numBlocks);
        simulationAvgVAFTheta = nan(1,numBlocks);
        simulationAvgVAFOmega = nan(1,numBlocks);
        simulationAvgVAFAlpha = nan(1,numBlocks);

        subjSDPos = nan(1,numBlocks);
        subjSDAcc = nan(1,numBlocks);
        subjSDTheta = nan(1,numBlocks);
        subjSDOmega = nan(1,numBlocks);
        subjSDAlpha = nan(1,numBlocks);

        subjAvgSDPos = nan(1,numBlocks);
        subjAvgSDAcc = nan(1,numBlocks);
        subjAvgSDTheta = nan(1,numBlocks);
        subjAvgSDOmega = nan(1,numBlocks);
        subjAvgSDAlpha = nan(1,numBlocks);

        simulationSDPos = nan(1,numBlocks);
        simulationSDAcc = nan(1,numBlocks);
        simulationSDTheta = nan(1,numBlocks);
        simulationSDOmega = nan(1,numBlocks);
        simulationSDAlpha = nan(1,numBlocks);
    end
    
    % Loop through blocks
    for blockNum = blockStart:blockEnd

        % Get information about experimental block
        [subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);
        
        % Load trimming indices
        load(blockStr,'Expression1')
        
        % Determine if current block is block 3 or block 4
        block3or4 = blockNumFunc(blockNum);

        % Folder holding best-fit trials from current block
        if(isfolder(strcat("best fit simulation/",parentFolder,"S10_B3")))
            % Folders named by subject number & block number
            simFolder = strcat("best fit simulation/",parentFolder,subjNum,"_B",block3or4);
        else
            % Folders named based on subject initials & trial date (old convention):
            simFolder = strcat("best fit simulation/",parentFolder,subjStr,trialDate,trialStr);
            simFolder = erase(simFolder, "_trial_");
        end
        

        % Folder holding trimmed experimental trials and averages
        expFolder = strcat("experimental data trimmed/",subjStr,trialDate,trialStr);
        expFolder = erase(expFolder, "_trial_");
        
        % Loop through trials, plotting each one and storing values in array
        for num = numStart:1:numEnd
            
            % Skip invalid trials
            if ismember(num,invalidTrials)
                VAFVel(num,blockNum) = NaN;
                
                % DEBUG
                disp(['Invalid trial:',num2str(num)])
                
                if ~justVel
                    VAFPos(num,blockNum) = NaN;
                    VAFAcc(num,blockNum) = NaN;
                    VAFTheta(num,blockNum) = NaN;
                    VAFOmega(num,blockNum) = NaN;
                    VAFAlpha(num,blockNum) = NaN;
                end
                continue
            end

            % Construct file name and load one experimental trial. Note
            % that experimental trials were trimmed and averages computed
            % using the script "statisticsExpMean.m"
            numStr = num2str(num);
            fileStr = strcat(expFolder,"/",numStr);
            load(fileStr,'pos','theta','vel','omega','acc','alpha','tdes',...
                'posAvg','thetaAvg','velAvg','omegaAvg','accAvg','alphaAvg');

            % Load best-fit simulated trial
            simProfiles = readmatrix(strcat(simFolder, "/best-fit profiles/", numStr, ".xlsx"));
            tdesSim     = simProfiles(end,1);
            posSim      = simProfiles(:,2);
            thetaSim    = simProfiles(:,3);
            velSim      = simProfiles(:,4);
            omegaSim    = simProfiles(:,5);
            accSim      = simProfiles(:,6);
            alphaSim    = simProfiles(:,7);
            
            % Trim experimental and simulated trajectories to same
            % duration, using shorter duration of the two
            if tdes > tdesSim
                len = int16(tdesSim/st+1);
                pos = pos(1:len);
                vel = vel(1:len);
                acc = acc(1:len);
                theta = theta(1:len);
                omega = omega(1:len);
                alpha = alpha(1:len);
                tdes = tdesSim;
            elseif tdesSim > tdes
                len = int16(tdes/st+1);
                posSim = posSim(1:len);
                velSim = velSim(1:len);
                accSim = accSim(1:len);
                thetaSim = thetaSim(1:len);
                omegaSim = omegaSim(1:len);
                alphaSim = alphaSim(1:len);
                tdesSim = tdes;
            end

            % Normalize experimental trajectories
            [pos,theta,vel,omega,acc,alpha] = normalize(pos,theta,vel,...
                omega,acc,alpha,tdes,st,stNorm,[]);
            
            % Normalize simulated trajectories
            [posSim,thetaSim,velSim,omegaSim,accSim,alphaSim] = ...
                normalize(posSim,thetaSim,velSim,omegaSim,accSim,alphaSim,...
                tdesSim,st,stNorm,[]);
            
            % Calculate simulated variable means
            velSimAvg      = mean(velSim);
            if ~justVel
                posSimAvg      = mean(posSim);
                thetaSimAvg    = mean(thetaSim);
                omegaSimAvg    = mean(omegaSim);
                accSimAvg      = mean(accSim);
                alphaSimAvg    = mean(alphaSim);
            end
            
            %%%% Compute Variance Accounted For (VAF) value for trial %%%%
            % Numerator sum of squares
            ssNumVel = sum((vel - velSim - velAvg + velSimAvg).^2);
            if ~justVel
                ssNumPos = sum((pos - posSim - posAvg + posSimAvg).^2);
                ssNumAcc = sum((acc - accSim - accAvg + accSimAvg).^2);
                ssNumTheta = sum((theta - thetaSim - thetaAvg + thetaSimAvg).^2);
                ssNumOmega = sum((omega - omegaSim - omegaAvg + omegaSimAvg).^2);
                ssNumAlpha = sum((alpha - alphaSim - alphaAvg + alphaSimAvg).^2);
            end
            
            % Denominator sum of squares
            ssDenVel = sum((vel - velAvg).^2);
            if ~justVel
                ssDenPos = sum((pos - posAvg).^2);
                ssDenAcc = sum((acc - accAvg).^2);
                ssDenTheta = sum((theta - thetaAvg).^2);
                ssDenOmega = sum((omega - omegaAvg).^2);
                ssDenAlpha = sum((alpha - alphaAvg).^2);
            end
            
            % Calculate VAF and store in array
            VAFVel(num,blockNum) = max(1 - ssNumVel/ssDenVel,0);
            if ~justVel
                VAFPos(num,blockNum) = max(1 - ssNumPos/ssDenPos,0);
                VAFAcc(num,blockNum) = max(1 - ssNumAcc/ssDenAcc,0);
                VAFTheta(num,blockNum) = max(1 - ssNumTheta/ssDenTheta,0);
                VAFOmega(num,blockNum) = max(1 - ssNumOmega/ssDenOmega,0);
                VAFAlpha(num,blockNum) = max(1 - ssNumAlpha/ssDenAlpha,0);
            end

        end
        
        % Calculate mean VAF value for one block
        blockAvgVAFVel(blockNum) = mean(VAFVel(:,blockNum),'omitnan');
        if ~justVel
            blockAvgVAFPos(blockNum) = mean(VAFPos(:,blockNum),'omitnan');
            blockAvgVAFAcc(blockNum) = mean(VAFAcc(:,blockNum),'omitnan');
            blockAvgVAFTheta(blockNum) = mean(VAFTheta(:,blockNum),'omitnan');
            blockAvgVAFOmega(blockNum) = mean(VAFOmega(:,blockNum),'omitnan');
            blockAvgVAFAlpha(blockNum) = mean(VAFAlpha(:,blockNum),'omitnan');
        end
        
        % Calulate mean and median VAF values for one subject across both blocks
        if ~mod(blockNum,2)
            subjAvgVAFVel(blockNum) = mean(blockAvgVAFVel(blockNum-1:blockNum),'omitnan');
            meanVAF(simNum,blockNum/2) = subjAvgVAFVel(blockNum);
            medianVAF(simNum,blockNum/2) = median(VAFVel(:,blockNum-1:blockNum),'all','omitnan');
                
            if ~justVel
                subjAvgVAFPos(blockNum) = mean(blockAvgVAFPos(blockNum-1:blockNum),'omitnan');
                subjAvgVAFAcc(blockNum) = mean(blockAvgVAFAcc(blockNum-1:blockNum),'omitnan');
                subjAvgVAFTheta(blockNum) = mean(blockAvgVAFTheta(blockNum-1:blockNum),'omitnan');
                subjAvgVAFOmega(blockNum) = mean(blockAvgVAFOmega(blockNum-1:blockNum),'omitnan');
                subjAvgVAFAlpha(blockNum) = mean(blockAvgVAFAlpha(blockNum-1:blockNum),'omitnan');
            end
        else
            subjAvgVAFVel(blockNum) = NaN;
            if ~justVel
                subjAvgVAFPos(blockNum) = NaN;
                subjAvgVAFAcc(blockNum) = NaN;
                subjAvgVAFTheta(blockNum) = NaN;
                subjAvgVAFOmega(blockNum) = NaN;
                subjAvgVAFAlpha(blockNum) = NaN;
            end
        end
        
        % Calculate standard deviation for one subject across both blocks
        if ~mod(blockNum,2)
            subjSDVel(blockNum) = std(VAFVel(:,blockNum-1:blockNum),1,'all','omitnan');
            sdevVAF(simNum,blockNum/2) = subjSDVel(blockNum);
            if ~justVel
                subjSDPos(blockNum) = std(VAFPos(:,blockNum-1:blockNum),1,'all','omitnan');
                subjSDAcc(blockNum) = std(VAFAcc(:,blockNum-1:blockNum),1,'all','omitnan');
                subjSDTheta(blockNum) = std(VAFTheta(:,blockNum-1:blockNum),1,'all','omitnan');
                subjSDOmega(blockNum) = std(VAFOmega(:,blockNum-1:blockNum),1,'all','omitnan');
                subjSDAlpha(blockNum) = std(VAFAlpha(:,blockNum-1:blockNum),1,'all','omitnan');
            end
        else
            subjSDVel(blockNum) = NaN;
            if ~justVel
                subjSDPos(blockNum) = NaN;
                subjSDAcc(blockNum) = NaN;
                subjSDTheta(blockNum) = NaN;
                subjSDOmega(blockNum) = NaN;
                subjSDAlpha(blockNum) = NaN;
            end
        end  
    end
    % End looping through blocks
    
    % Calculate mean VAF for all trials run using current simulation
    simulationAvgVAFVel(1) = mean(VAFVel,'all','omitnan');
    if ~justVel
        simulationAvgVAFPos(1) = mean(VAFPos,'all','omitnan');
        simulationAvgVAFAcc(1) = mean(VAFAcc,'all','omitnan');
        simulationAvgVAFTheta(1) = mean(VAFTheta,'all','omitnan');
        simulationAvgVAFOmega(1) = mean(VAFOmega,'all','omitnan');
        simulationAvgVAFAlpha(1) = mean(VAFAlpha,'all','omitnan');
    end
    
    % Calculate standard deviation of subject mean VAFs
    sdevSubjMeanVel(simNum) = std(meanVAF(simNum, :));
    subjAvgSDVel(1) = mean(subjSDVel,'all','omitnan');
    if ~justVel
        subjAvgSDPos(1) = mean(subjSDPos,'all','omitnan');
        subjAvgSDAcc(1) = mean(subjSDAcc,'all','omitnan');
        subjAvgSDTheta(1) = mean(subjSDTheta,'all','omitnan');
        subjAvgSDOmega(1) = mean(subjSDOmega,'all','omitnan');
        subjAvgSDAlpha(1) = mean(subjSDAlpha,'all','omitnan');
    end
    
    % Calculate standard deviation for all trials run using current simulation
    simulationSDVel(1) = std(VAFVel,1,'all','omitnan');
    if ~justVel
        simulationSDPos(1) = std(VAFPos,1,'all','omitnan');
        simulationSDAcc(1) = std(VAFAcc,1,'all','omitnan');
        simulationSDTheta(1) = std(VAFTheta,1,'all','omitnan');
        simulationSDOmega(1) = std(VAFOmega,1,'all','omitnan');
        simulationSDAlpha(1) = std(VAFAlpha,1,'all','omitnan');
    end
    
    % Calculate mean and standard deviation of subject median VAFs
    meanOfMedianVAF(simNum) = mean(medianVAF(simNum,:));
    sdevMedianVAF(simNum) = std(medianVAF(simNum,:));
    
    % Add values to table saving mean & SD for all simulations
    overallResults(:,simNum-simStart+2) = {simLabels(simNum);...
        simulationAvgVAFVel(1);sdevSubjMeanVel(simNum);meanOfMedianVAF(simNum);...
        sdevMedianVAF(simNum)};
    
    if saveFile

        if(not(isfolder(strcat("best fit simulation/",parentFolder,"_Statistics"))))
            mkdir(strcat("best fit simulation/",parentFolder),"_Statistics");
        end
        saveFolder = strcat("best fit simulation/",parentFolder,"/_Statistics/");

        % Define column names
        blockHeaders = {'S1 Block 3','S1 Block 4','S2 Block 3',...
            'S2 Block 4','S3 Block 3','S3 Block 4','S4 Block 3',...
            'S4 Block 4','S5 Block 3','S5 Block 4','S6 Block 3',...
            'S6 Block 4','S7 Block 3','S7 Block 4','S8 Block 3',...
            'S8 Block 4','S9 Block 3','S9 Block 4','S10 Block 3',...
            'S10 Block 4','S11 Block 3','S11 Block 4'};
        blockHeaders = blockHeaders(blockStart:blockEnd);
        
        % Define row names
        rowNames = string(numStart:1:numEnd);
        
        % Put parameter arrays in cell array so for loop can select them
        if ~justVel
            VAFArray = {VAFPos, VAFVel, VAFAcc,...
                VAFTheta, VAFOmega, VAFAlpha};
            blockAvgVAFArray = {blockAvgVAFPos, blockAvgVAFVel, blockAvgVAFAcc,...
                blockAvgVAFTheta, blockAvgVAFOmega, blockAvgVAFAlpha};
            subjAvgVAFArray = {subjAvgVAFPos, subjAvgVAFVel, subjAvgVAFAcc,...
                subjAvgVAFTheta, subjAvgVAFOmega, subjAvgVAFAlpha};
            simulationAvgVAFArray = {simulationAvgVAFPos, simulationAvgVAFVel,...
                simulationAvgVAFAcc, simulationAvgVAFTheta,...
                simulationAvgVAFOmega, simulationAvgVAFAlpha};
            subjSDArray = {subjSDPos, subjSDVel, subjSDAcc, subjSDTheta,...
                subjSDOmega, subjSDAlpha};
            subjAvgSDArray = {subjAvgSDPos, subjAvgSDVel,subjAvgSDAcc,...
                 subjAvgSDTheta, subjAvgSDOmega,subjAvgSDAlpha};
            simulationSDArray = {simulationSDPos, simulationSDVel, simulationSDAcc,...
                simulationSDTheta, simulationSDOmega, simulationSDAlpha};
        
            paramNameArray = {'1. Position', '2. Velocity', '3. Acceleration',...
                '4. Angle', '5. Angular Vel', '6. Angular Acc'};
        
            for param = 1:6
                VAFTable = array2table(VAFArray{param},...
                    'VariableNames',blockHeaders,'RowNames',rowNames);    
                blockAvgVAFTable = array2table(blockAvgVAFArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Mean VAF per Block"]);
                subjAvgVAFTable = array2table(subjAvgVAFArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Mean VAF per Subject"]);
                simAvgVAFTable = array2table(simulationAvgVAFArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Mean VAF for Simulation"]);
                subjSDTable = array2table(subjSDArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Std Dev per Subject"]);
                subjAvgSDTable = array2table(subjAvgSDArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Mean of Per-Subject Std Devs"]);
                simulationSDTable = array2table(simulationSDArray{param},...
                    'VariableNames',blockHeaders,'RowNames',["Overall Std Dev for Simulation"]);

                outputTable = [VAFTable; blockAvgVAFTable; subjAvgVAFTable;...
                    simAvgVAFTable; subjSDTable; subjAvgSDTable; simulationSDTable];

                paramName = paramNameArray{param};
                parentFolder = erase(parentFolder, "/");
                fileName = strcat("best fit simulation/",parentFolder,"/_Statistics/", paramName, " R Squared - ",parentFolder,"2022-01-14.xlsx");
                writetable(outputTable,fileName,'WriteRowNames',true);    
            end
        end
        
    % Save results in excel table and .mat file
    resultsTable = cell2table(overallResults);
    writetable(resultsTable,"best fit simulation/overall VAFs.xlsx",'WriteVariableNames',false);
    save('best fit simulation/VAFarrays.mat','meanVAF','sdevVAF', 'medianVAF', 'sdevMedianVAF')

    end
    
end

%{
% Statistics
% original_medians = medianVAF(26,:);
multimode_medians = medianVAF(30,:);
rigid_body_medians = medianVAF(29,:);

% % Check if rigid-body model VAF is significantly higher than original input
% % shaping VAF
% [T, sig5, s25, s1] = sigDiffPaired(rigid_body_medians, original_medians)

% Check if rigid-body model VAF is significantly higher than multi-mode
% input shaping VAF
[T, sig5, s25, s1] = sigDiffPaired(rigid_body_medians, multimode_medians)

% % Check if multi-mode input shaping VAF is significantly is significantly
% % higher than original input shaping VAF
% [T, sig5, s25, s1] = sigDiffPaired(multimode_medians, original_medians)
%}

toc