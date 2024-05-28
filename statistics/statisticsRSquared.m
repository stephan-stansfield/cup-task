% statisticsRSquared
    
clear all;
close all;
clc;

addpath('data', 'trim times', genpath('experimental data'));
st          = 0.001;                                                        % Simulation time step (s)
stNorm      = 0.0001;                                                       % Normalized time step (s)
tNorm       = 0:stNorm:1;                                                   % Normalized time vector
normLen     = length(tNorm);    

% Choose whether to take average of all experimental and simulated
% trials in block, or to take early and late simulated trials only.
type = "all";
%     type = "early and late";

% Choose whether to save generated statistics
save = true;

% Choose which experimental blocks to visualize (all acceptable blocks
% would be 1-22)
blockStart  = 1;
blockEnd    = 22;
numBlocks   = blockEnd - blockStart + 1;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Cycle through simulation types
for count = 1
    % 1. Nominal Input Shaping (No Impedance Model)
    % 2. Multi-Mode Input Shaping, FF
    % 3. Multi-Mode Input Shaping, No FF
    % 4. Slow Mode IS, FF
    % 5. Slow Mode IS, No FF
    % 6. No Impedance IS, FF
    % 7. No Impedance IS, No FF
    % 8. Fast Mode IS, FF
    % 9. Fast Mode IS, No FF
    % 10. Rigid Body IS, FF
    % 11. Rigid Body IS, No FF
    % 12. Slow Mode IS, FF, Velocity Weights Only
    % 13. Fast Mode IS, FF, Velocity Weights Only
    
    % Get name of folder holding best-fit simulations for selected
    % simulation type
    parentFolder = simulationFolder(count);
    if parentFolder == "/"
        continue
    end
    
    % Initialize arrays to hold r squared values
    rSquaredPos = zeros(50,numBlocks);
    rSquaredVel = zeros(50,numBlocks);
    rSquaredAcc = zeros(50,numBlocks);
    rSquaredTheta = zeros(50,numBlocks);
    rSquaredOmega = zeros(50,numBlocks);
    rSquaredAlpha = zeros(50,numBlocks);
%     ssTotArray = zeros(1,50);
%     ssResArray = zeros(1,50);

    blockAvgRSqPos = zeros(1,numBlocks);
    blockAvgRSqVel = zeros(1,numBlocks);
    blockAvgRSqAcc = zeros(1,numBlocks);
    blockAvgRSqTheta = zeros(1,numBlocks);
    blockAvgRSqOmega = zeros(1,numBlocks);
    blockAvgRSqAlpha = zeros(1,numBlocks);
    
    subjAvgRSqPos = zeros(1,numBlocks);
    subjAvgRSqVel = zeros(1,numBlocks);
    subjAvgRSqAcc = zeros(1,numBlocks);
    subjAvgRSqTheta = zeros(1,numBlocks);
    subjAvgRSqOmega = zeros(1,numBlocks);
    subjAvgRSqAlpha = zeros(1,numBlocks);
    
    simulationAvgRSqPos = zeros(1,numBlocks);
    simulationAvgRSqVel = zeros(1,numBlocks);
    simulationAvgRSqAcc = zeros(1,numBlocks);
    simulationAvgRSqTheta = zeros(1,numBlocks);
    simulationAvgRSqOmega = zeros(1,numBlocks);
    simulationAvgRSqAlpha = zeros(1,numBlocks);
    
%     if ~mod(numBlocks,2)
%         subjAvgRSq = zeros(1,numBlocks/2);
%     end
    
    for blockNum = blockStart:blockEnd

        % Get information about experimental block
        [subjNum, subjStr, trialDate, trialStr, blockStr, ~,...
            minVelOutliers, invalidTrials] = blockDictionary(blockNum);
        
        % Load trimming indices
        load(blockStr,'Expression1')

        % Folder holding best-fit trials from current block
        simFolder = strcat("best fit simulation/",parentFolder,subjStr,trialDate,trialStr);
        simFolder = erase(simFolder, "_trial_");

        % Folder holding trimmed experimental trials and averages
        expFolder = strcat("experimental data trimmed/",subjStr,trialDate,trialStr);
        expFolder = erase(expFolder, "_trial_");
        
        % Loop through trials, plotting each one and storing values in array
        for num = numStart:1:numEnd
            
            % Skip invalid trials
            if ismember(num,invalidTrials)
                rSquaredPos(num,blockNum) = NaN;
                rSquaredVel(num,blockNum) = NaN;
                rSquaredAcc(num,blockNum) = NaN;
                rSquaredTheta(num,blockNum) = NaN;
                rSquaredOmega(num,blockNum) = NaN;
                rSquaredAlpha(num,blockNum) = NaN;
                continue
            end

            % Construct file name and load one experimental trial
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
            
             % Trim experimental and simulated trajectories to same duration
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
            [pos,theta,vel,omega,acc,alpha] = normalize(pos,theta,vel,omega,...
                acc,alpha,tdes,st,stNorm,[]);
            
            % Normalize simulated trajectories
            [posSim,thetaSim,velSim,omegaSim,accSim,alphaSim] = normalize(posSim,...
                thetaSim,velSim,omegaSim,accSim,alphaSim,tdesSim,st,stNorm,[]);
            
            %%%% Compute Coefficient of Determination (r^2) value for trial
            % Total sum of squares
            ssTotPos = sum((pos - posAvg).^2);
            ssTotVel = sum((vel - velAvg).^2);
            ssTotAcc = sum((acc - accAvg).^2);
            ssTotTheta = sum((theta - thetaAvg).^2);
            ssTotOmega = sum((omega - omegaAvg).^2);
            ssTotAlpha = sum((alpha - alphaAvg).^2);
            
            % Residual sum of squares
            ssResPos = sum((pos - posSim).^2);
            ssResVel = sum((vel - velSim).^2);
            ssResAcc = sum((acc - accSim).^2);
            ssResTheta = sum((theta - thetaSim).^2);
            ssResOmega = sum((omega - omegaSim).^2);
            ssResAlpha = sum((alpha - alphaSim).^2);
            
            % Calculate simple coefficient of determination and store in
            % array
            rSquaredPos(num,blockNum) = 1 - ssResPos/ssTotPos;
            rSquaredVel(num,blockNum) = 1 - ssResVel/ssTotVel;
            rSquaredAcc(num,blockNum) = 1 - ssResAcc/ssTotAcc;
            rSquaredTheta(num,blockNum) = 1 - ssResTheta/ssTotTheta;
            rSquaredOmega(num,blockNum) = 1 - ssResOmega/ssTotOmega;
            rSquaredAlpha(num,blockNum) = 1 - ssResAlpha/ssTotAlpha;
            
%             % Regression sum of squares
%             ssRegPos = sum((posSim-posAvg).^2);
%             ssRegVel = sum((velSim-velAvg).^2);
%             ssRegAcc = sum((accSim-accAvg).^2);
%             ssRegTheta = sum((thetaSim-thetaAvg).^2);
%             ssRegOmega = sum((omegaSim-omegaAvg).^2);
%             ssRegAlpha = sum((alphaSim-alphaAvg).^2);
            
%             % Store SS values in arrays to examine them
%             ssTotArray(num) = ssTotVel;
%             ssResArray(num) = ssResVel;

        end
        
        % Take average r^2 value for one block
        blockAvgRSqPos(blockNum) = mean(rSquaredPos(:,blockNum),'omitnan');
        blockAvgRSqVel(blockNum) = mean(rSquaredVel(:,blockNum),'omitnan');
        blockAvgRSqAcc(blockNum) = mean(rSquaredAcc(:,blockNum),'omitnan');
        blockAvgRSqTheta(blockNum) = mean(rSquaredTheta(:,blockNum),'omitnan');
        blockAvgRSqOmega(blockNum) = mean(rSquaredOmega(:,blockNum),'omitnan');
        blockAvgRSqAlpha(blockNum) = mean(rSquaredAlpha(:,blockNum),'omitnan');
        
        % Take average r^2 value for one subject across both blocks
        if ~mod(blockNum,2)
            subjAvgRSqPos(blockNum) = mean(blockAvgRSqPos(blockNum-1:blockNum),'omitnan');
            subjAvgRSqVel(blockNum) = mean(blockAvgRSqVel(blockNum-1:blockNum),'omitnan');
            subjAvgRSqAcc(blockNum) = mean(blockAvgRSqAcc(blockNum-1:blockNum),'omitnan');
            subjAvgRSqTheta(blockNum) = mean(blockAvgRSqTheta(blockNum-1:blockNum),'omitnan');
            subjAvgRSqOmega(blockNum) = mean(blockAvgRSqOmega(blockNum-1:blockNum),'omitnan');
            subjAvgRSqAlpha(blockNum) = mean(blockAvgRSqAlpha(blockNum-1:blockNum),'omitnan');
        else
            subjAvgRSqPos(blockNum) = 0;
            subjAvgRSqVel(blockNum) = 0;
            subjAvgRSqAcc(blockNum) = 0;
            subjAvgRSqTheta(blockNum) = 0;
            subjAvgRSqOmega(blockNum) = 0;
            subjAvgRSqAlpha(blockNum) = 0;
        end
        
    end
    
    % Calculate average r^2 for all trials run using current simulation
    simulationAvgRSqPos(1) = mean(rSquaredPos,'all','omitnan');
    simulationAvgRSqVel(1) = mean(rSquaredVel,'all','omitnan');
    simulationAvgRSqAcc(1) = mean(rSquaredAcc,'all','omitnan');
    simulationAvgRSqTheta(1) = mean(rSquaredTheta,'all','omitnan');
    simulationAvgRSqOmega(1) = mean(rSquaredOmega,'all','omitnan');
    simulationAvgRSqAlpha(1) = mean(rSquaredAlpha,'all','omitnan');
    
    if save

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
        rSquaredArray = {rSquaredPos, rSquaredVel, rSquaredAcc,...
            rSquaredTheta, rSquaredOmega, rSquaredAlpha};
        blockAvgRSqArray = {blockAvgRSqPos, blockAvgRSqVel, blockAvgRSqAcc,...
            blockAvgRSqTheta, blockAvgRSqOmega, blockAvgRSqAlpha};
        subjAvgRSqArray = {subjAvgRSqPos, subjAvgRSqVel, subjAvgRSqAcc,...
            subjAvgRSqTheta, subjAvgRSqOmega, subjAvgRSqAlpha};
        simulationAvgRSqArray = {simulationAvgRSqPos, simulationAvgRSqVel,...
            simulationAvgRSqAcc, simulationAvgRSqTheta,...
            simulationAvgRSqOmega, simulationAvgRSqAlpha};
        paramNameArray = {'1. Position', '2. Velocity', '3. Acceleration',...
            '4. Angle', '5. Angular Vel', '6. Angular Acc'};
        
        
        for param = 1:6
        
            rSqTable = array2table(rSquaredArray{param},...
            'VariableNames',blockHeaders,'RowNames',rowNames);    
            blockAvgSqTable = array2table(blockAvgRSqArray{param},...
                'VariableNames',blockHeaders,'RowNames',["Average R^2 for Block"]);
            subjAvgSqTable = array2table(subjAvgRSqArray{param},...
                'VariableNames',blockHeaders,'RowNames',["Average R^2 for Subject"]);
            simAvgSqTable = array2table(simulationAvgRSqArray{param},...
                'VariableNames',blockHeaders,'RowNames',["Overall R^2 for Simulation"]);

            outputTable = [rSqTable; blockAvgSqTable; subjAvgSqTable; simAvgSqTable];
            
%             rSqTable = array2table(rSquaredVel,'VariableNames',blockHeaders,...
%                 'RowNames',rowNames);    
%             blockAvgSqTable = array2table(blockAvgRSq,'VariableNames',blockHeaders,...
%                 'RowNames',["Average R^2 for Block"]);
%             subjAvgSqTable = array2table(subjAvgRSq,'VariableNames',blockHeaders,...
%                 'RowNames',["Average R^2 for Subject"]);
%             simAvgSqTable = array2table(simulationAvgRSq,'VariableNames',blockHeaders,...
%                 'RowNames',["Overall R^2 for Simulation"]);
% 
%             outputTable = [rSqTable; blockAvgSqTable; subjAvgSqTable; simAvgSqTable];

            paramName = paramNameArray{param};
            parentFolder = erase(parentFolder, "/");
            fileName = strcat("best fit simulation/",parentFolder,"/_Statistics/", paramName, " R Squared - ",parentFolder,".xlsx");
            writetable(outputTable,fileName,'WriteRowNames',true);    
        
        end
    end
    
end
        