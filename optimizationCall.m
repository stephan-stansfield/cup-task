% OPTIMIZATIONCALL
%
% This script allows choice of parameters and of blocks/trials to input
% to the optimization.m function. That function runs an optimization to
% find the simulated trials that best fit the experimental data .

close all
clear

addpath('data', genpath('experimental data'));

% Max number of evaluations per trial optimization algorithm will compute
maxEval  = 10;

% Choose range of experimental trials to fit (trial numbers start at 0)
numStart    = 0;
numEnd      = 49;

% Choose blocks to loop through. The "testing blocks" 3 and 4 for all 11 
% subjects correspond to blocks 1-22
blockStart  = 1;
blockEnd    = 2;

% Set range of additional duration (in seconds) to add to simulated trials
% before trimming to start and stop thresholds
delayMin = 0;
delayMax = 1.000;

% Choose whether to plot individual trial velocities (true) or to plot and
% save whole blocks of trials (false, the default)
plotInd = false;

% Plot desired velocity trajectory (only applied when printing an
% individual trial)
printDes = false;

% Print system transfer function, poles, and zeros
printSys = false;

% Plot simulated trajectories, the default (set this to false to only 
% plot experimental trajectories)
plotSim = true;

% Choose whether to identify location of velocity peaks and valleys on
% plots
plotPeaks = false;

% Choose whether to fit optimization parameters to each individual
% trial (true, the default), or to fit one set of parameters to an
% entire block
fitEachTrial = true;

% Choose optimization type
optimizationType = "input shaping 2 impulse impedance";
%     optimizationType = "input shaping 2 impulse no impedance";
%     optimizationType = "input shaping 4 impulse";

% Designate type of internal model used for input shaping and feedforward
%     intModel = "full";
    intModel = "rigid body";
%     intModel = "no impedance";
%     intModel = "slow";
%     intModel = "fast";

% Choose if simulated (external) model includes hand impedance
impedance = false;

% Choose linear or nonlinear simulation
simVersion  = "nonlinear";
% simVersion  = "linear";

% Set kinematic variable weights for RMSE objective function. Order of
% array is pos, theta, vel, omega, acc, alpha
weights = [10, 10, 5, 5, 1, 1];

% Designate model with or without feedforward force term
forwardF = true;

% If running multiple iterations in a row, can use counter and if
% statements to designate desired parameters
for setting = 11
    
    if setting == 1
        % 1. Nominal input shaping without hand impedance
        optimizationType = "input shaping 2 impulse no impedance";
        intModel = "no impedance";
        impedance = false;
        simVersion = "nonlinear";
        forwardF = true;

    elseif setting == 2
        % 2. Multi-mode internal model, no feedforward force
        optimizationType = "input shaping 4 impulse";
        intModel = "full";
        impedance = true;
        simVersion  = "nonlinear";
        forwardF = false;

    elseif setting == 3
        % 3. Slow mode internal model, no feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "slow";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = false;

    elseif setting == 4
        % 4. Fast mode internal model, no feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "fast";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = false;

    elseif setting == 5
        % 5. Rigid body internal model, no feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "rigid body";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = false;
        
    elseif setting == 6
        % 6. No impedance internal model, no feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "no impedance";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = false;
        
    elseif setting == 7
        % 7. Multi-Mode internal model, feedforward force
        optimizationType = "input shaping 4 impulse";
        intModel = "full";
        impedance   = true;
        simVersion = "nonlinear";
        forwardF = true;

    elseif setting == 8
        % 8. Slow mode internal model, feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "slow";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = true;
        
    elseif setting == 9
        % 9. Fast mode internal model, feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "fast";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = true; 

        elseif setting == 10
        % 10. Rigid body internal model, feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "rigid body";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = true;          
        
    elseif setting == 11
        % 11. No impedance internal model, feedforward force
        optimizationType = "input shaping 2 impulse impedance";
        intModel = "no impedance";
        impedance = true;
        simVersion = "nonlinear";
        forwardF = true;
        
    
    end
            
    % Start timing total code execution time
    tic

    % Generate name of subfolder where data will be stored
    if optimizationType == "input shaping 4 impulse"
        oString = '4 Impulse, ';
        rangeEnd = "P";
    elseif optimizationType == "input shaping 2 impulse impedance"
        oString = '2 Impulse, Impedance, ';
        rangeEnd = "L";
    elseif optimizationType == "input shaping 2 impulse no impedance"
        oString = '2 Impulse, No Impedance, ';
        rangeEnd = "J";
    end

    if forwardF
        fString = 'Feedforward F';
    else
        fString = 'No Feedforward F';
    end
    
    if intModel == "rigid body"
        mString = ', Rigid Body I.S. Simplification';
    elseif intModel == "no impedance"
        mString = ', No Impedance I.S. Simplification';
    elseif intModel == "slow"
        mString = ', Slow Mode I.S. Simplification';
    elseif intModel == "fast"
        mString = ', Fast Mode I.S. Simplification';
    else
        mString = '';
    end
        
    modelType = [oString, fString, mString];

    nameStr = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
    modelStr = strcat(nameStr,"_",modelType);
    mkdir("data", modelStr);

    % Loop through selected blocks
    for blockNum = blockStart:blockEnd
        
        optimization(blockNum, optimizationType, modelStr,forwardF,...
            intModel, impedance, simVersion, weights, delayMin,...
            delayMax, plotInd, printDes, printSys, plotSim, plotPeaks,...
            maxEval, numStart, numEnd)

        % Close all windows at the end of running a block (necessary if
        % running multiple blocks)
        if ~plotInd || ~printDes
            close all
        end

    end

    % Name folder that will hold all data for this iteration
    topFolderStr = strcat("data/", modelStr, "/");

    % Write data from all blocks into a single spreadsheet
    filename = strcat(topFolderStr,modelType,".xlsx");
    
    if optimizationType == "input shaping 4 impulse" ||...
            optimizationType == "input shaping 2 impulse impedance"
        variableNames = {'Trial','Fmin','VRMSE','B','K','Tdelay'};
    elseif optimizationType == "input shaping 2 impulse no impedance"
        variableNames = {'Trial','Fmin','VRMSE','Tdelay'};
    end
    
    % Write headers to spreadsheet
    writecell(variableNames,filename,'Sheet',1)
    
    % Load data for each block and write to spreadsheet
    for  blockNum = blockStart:1:blockEnd
        load(strcat(topFolderStr, num2str(blockNum), ".mat"),"returnArray");
        if plotInd
            numStart = 1;
            numEnd = 1;
        else
            excelRowStart = (blockNum-1)*50 + 2;
            excelRowEnd = excelRowStart + 49;
        end
        range = strcat("A",num2str(excelRowStart),":",rangeEnd,num2str(excelRowEnd));
        writematrix(returnArray,filename,'Sheet',1,'Range',range)
    end

    % Stop timing code execution
    toc
    
end