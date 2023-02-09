% This function allows choice of parameters and of blocks/trials to input
% to the optimization.m function. That function runs an optimization to
% find the simulated trials that best fit the experimental data .

close all
clear all

addpath('data', genpath('experimental data'));

% Max number of evaluations per trial optimization algorithm will compute
setMaxEval  = 100;

% Choose range of experimental trials to fit (trial numbers start at 0)
numStart    = 0;
numEnd      = 0;

% Choose blocks to loop through. The 11 acceptable subjects 
% correspond to blocks 1-22. All blocks would be 1 to 28.
blockStart  = 1;
blockEnd    = 1;

% Set range of additional duration values to simulate trials before trimming
delayMin = 0;       % [sec]
delayMax = 0.250;   % [sec]

% Choose whether to plot and save whole blocks of trials (false -- 
% the default) or to plot individual trial velocities (true)
plotInd = false;

% Plot desired velocity trajectory (only applied when printing an
% individual trial)
printDes = false;

% Print system transfer function, poles, and zeros
printSys = false;

% Plot simulated trajectories, the default (set this to false to *only* 
% plot experimental trajectories)
plotSim = true;

% Choose whether to identify location of velocity peaks and valleys on
% plots
plotPeaks = false;

% Choose whether to fit optimization parameters to each individual
% trial (true - the default), or to fit one set of parameters to an
% entire block
fitEachTrial = true;

% Choose optimization type
optimizationType = "input shaping 2 impulse impedance";
%     optimizationType = "input shaping 2 impulse no impedance";
%     optimizationType = "input shaping 4 impulse";
%     optimizationType = "submovement";

% Designate type of internal model used for input shaping and feedforward.
%     ver = "full";
%     ver = "rigid body";
    ver = "no impedance";
%     ver = "slow";
%     ver = "fast";

% Choose if simulation includes time and/or amplitude error
timeMod = false;
ampMod  = false;

% Choose if simulated model includes hand impedance
impedance = true;

% Choose linear, nonlinear, or linear + pendulum lock simulation
simVersion  = "pendLock";
%     simVersion  = "linear";
%     simVersion  = "nonlinear";

% Designate objective function type
objective = "RMSE";
%     objective = "features";

% Set kinematic variable weights for RMSE objective function. Order of
% array is pos, theta, vel, omega, acc, alpha (in future, make this a
% structure)
weights = [10, 10, 5, 5, 1, 1];
%     weights = [0, 0, 1, 0, 0, 0];

% Designate model with feedforward F or with no feedforward F term
forwardF = true;

% Choose whether to fit optimization parameters to each individual
% trial or to fit to the average of an entire block of trials
fitMethod = "eachTrial";
%     fitMethod = "fitToAverage";

% If running multiple iterations in a row, can use counter and if
% statements to designate desired parameters.
for setting = 0
    
    if setting == 1
        % 1. Original input shaping
        optimizationType = "input shaping 2 impulse no impedance";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = false;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;

    elseif setting == 2
        % 2. Multi-mode, no FF
        optimizationType = "input shaping 4 impulse";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;

    elseif setting == 3
        % 3. Slow mode, no FF
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "slow";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;

    elseif setting == 4
        % 4. Fast mode, no FF
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "fast";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;

    elseif setting == 5
        % 5. Rigid body, no FF
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "rigid body";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;
        
    elseif setting == 6
        % 6. No impedance, no FF
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "no impedance";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;
        
    elseif setting == 7
        % 7. Multi-Mode, Feedforward
        optimizationType = "input shaping 4 impulse";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;

    elseif setting == 8
        % 8. Slow mode, Feedforward
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "slow";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;
        
    elseif setting == 9
        % 9. Fast mode, feedforward
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "fast";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true; 
        
    elseif setting == 10
        % 10. Rigid body simplification, feedforward
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "rigid body";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;          
        
    elseif setting == 11
        % 11. No impedance simplification, feedforward
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "no impedance";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif setting == 12
        % 12. Original input shaping, feature-based objective
        optimizationType = "input shaping 2 impulse no impedance";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = false;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;
        
    elseif setting == 13
        % 13. Multi-Mode, Feedforward, feature-based objective
        optimizationType = "input shaping 4 impulse";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;
          
    elseif setting == 14
        % 14. No impedance simplification, feedforward, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "no impedance";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;
        
    elseif setting == 15
        % 15. Rigid body simplification, feedforward, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "rigid body";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;
        
    elseif setting == 16
        % 16. Slow mode, Feedforward, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "slow";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;
        
    elseif setting == 17
        % 17. Fast mode, feedforward, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "fast";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = true;            
        
    elseif setting == 18
        % 18. Multi-mode, no FF, feature-based objective
        optimizationType = "input shaping 4 impulse";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = false;
        
    elseif setting == 19
        % 19. No impedance, no FF, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "no impedance";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = false;
        
    elseif setting == 20
        % 20. Rigid body, no FF, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "rigid body";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = false;
        
    elseif setting == 21
        % 21. Slow mode, no FF, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "slow";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = false;
        
    elseif setting == 22
        % 22. Fast mode, no FF, feature-based objective
        optimizationType = "input shaping 2 impulse impedance";
        ver         = "fast";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "features";
        forwardF = false;
        
    elseif setting == 23
        % 23. Submovement Input, No Impedance
        optimizationType = "submovement";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = false;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;
        
    elseif setting == 24
        % 24. Submovement Input, Impedance, Forward F
        optimizationType = "submovement";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = true;
        
    elseif setting == 25
        % 25. Submovement Input, Impedance, No Forward F
        optimizationType = "submovement";
        ver         = "";
        timeMod     = false;
        ampMod      = false;
        impedance   = true;
        simVersion  = "pendLock";
        objective = "RMSE";
        forwardF = false;
        
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
    elseif optimizationType == "submovement"
        oString = 'Submovements, ';
        rangeEnd = "K";
    end

    if timeMod
        tString = 'Var Timing, ';
    else
        tString = 'No Var Timing, ';
    end

    if ampMod
        aString = 'Var Amp, ';
    else
        aString = 'No Var Amp, ';
    end

    if forwardF
        fString = 'Feedforward F';
    else
        fString = 'No Feedforward F';
    end
    
    if ver == "rigid body"
        mString = ', Rigid Body I.S. Simplification';
    elseif ver == "no impedance"
        mString = ', No Impedance I.S. Simplification';
    elseif ver == "slow"
        mString = ', Slow Mode I.S. Simplification';
    elseif ver == "fast"
        mString = ', Fast Mode I.S. Simplification';
    else
        mString = '';
    end
        
    modelType = [oString, tString, aString, fString, mString];

    nameStr     = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
    modelStr    = strcat(nameStr,"_",modelType);
    mkdir("data", modelStr);
    
    if fitMethod == "fitToAverage"
        plotInd = true;
    end

    % Loop through selected blocks
    for blockNum = blockStart:blockEnd
        
        % Selects type of fitting (per trial versus per block)
        if fitEachTrial
            optimization(blockNum,optimizationType,modelStr,...
                forwardF,timeMod,ampMod,ver,impedance,simVersion,...
                objective,weights,delayMin,delayMax,plotInd,printDes,...
                printSys,plotSim,plotPeaks,setMaxEval,numStart,...
                numEnd,fitMethod)

        else
            optimization_block_20210106(subjStr,trialDate,trialStr,blockStr,...
                blockNum,optimizationType,modelStr,forwardF,timeMod,ampMod,ver)
        end

        % Close all windows at the end of running a block (practically
        % necessary if running multiple blocks)
        if ~plotInd || ~printDes
%             if ~printDes
            close all
%             end
        end

    end

    % Name folder that will hold all data for this iteration
    topFolderStr = strcat("data/", modelStr, "/");

    % Write data from all blocks into a single spreadsheet
    
%         if ~plotInd
        filename = strcat(topFolderStr,modelType,".xlsx");
        
        if optimizationType == "input shaping 4 impulse"
            variableNames = {'Trial','Fmin','VRMSE','B','K','Impulse1','Impulse2',...
                            'Impulse3','Impulse4','Tdelay','Amp11','Amp12','Amp21',...
                            'Amp22','Sub1End','Sub2Start'};
        elseif optimizationType == "input shaping 2 impulse impedance"
            variableNames = {'Trial','Fmin','VRMSE','B','K','Impulse1','Impulse2',...
                            'Tdelay','Amp1','Amp2','Sub1End','Sub2Start'};
        elseif optimizationType == "input shaping 2 impulse no impedance"
            variableNames = {'Trial','Fmin','VRMSE','Impulse1','Impulse2',...
                            'Tdelay','Amp1','Amp2','Sub1End','Sub2Start'};
        elseif optimizationType == "submovement"
            variableNames = {'Trial','Fmin','VRMSE','B','K','Dist1',...
                            'Submovement 1 end','Submovement 2 start','Tdelay'};
        end
        
        % Write headers to spreadsheet
        writecell(variableNames,filename,'Sheet',1)
        
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
%         end

    % Stop timing code execution
    toc
    
end

% Play a sound to indicate code is finished running
%     load handel
%     sound(y,Fs)
beep