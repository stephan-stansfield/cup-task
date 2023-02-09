% optimization
%
% This function uses the NLopt plugin to run a nonlinear optimization
% algorithm and fit simulated cup task motion profiles to experimentally
% measured human motion profiles. Simulated motions are generated using a
% selected combination of:
%
% * Input profile generation (external function call: simInputShape).
%   Either:
%   ** Multi-mode input shaping (4 impulses)
%   ** Single-mode input shaping
%       *** 2 impulses based on rigid-body simplification (lumped mass)
%       *** 2 impulses based on no-impedance simplification (just pendulum)
%       *** 2 impulses based on slower mode of full model
%       *** 2 impulses based on faster mode of full model
%   ** Two pure minimum jerk submovements
%
% * Physical plant
%   ** Including vs. not including hand impedance
%   ** Linear approximation vs. nonlinear model
%
% * Control structure
%   ** Including vs. not including feedforward force term
%
% Note that NLopt(args) is a nested function within this function. NLopt
% sets the optimization parameters and calls the optimization algorithm.
% objFunc.m is the objective function which NLopt(args) calls repeatedly.
%
% This function calls external functions sysCreate.m, simInputShape.m, and
% trimData.m. It also loads external files containing the experimental
% profiles and proper indices at which to trim the data. These files must
% be in the MATLAB path in order for this function to run, which the
% function optimizationCall.m should automatically do. (Adds folders "data"
% and "trim times" to the path, as well as "experimental data" and its
% subfolders.)

function optimization(blockNum,optimizationType,modelStr,forwardF,timeMod,...
        ampMod,ver,impedance,simVersion,objective, weights,delayMin,delayMax,...
        plotInd,printDes,printSys,plotSim,plotPeaks,setMaxEval,numStart,...
        numEnd,fitMethod)

% Set variable values (note that constants for gravity, physical system 
% mass, and length are defined in globalData.m)
Bmin        = 7.5;                                                            % Minimum damping value (N*s/m)
Bmax        = 75;                                                           % Maximum damping value (N*s/m)
Kmin        = 7.5;                                                            % Minimum  stiffness value (N/m)
Kmax        = 750;                                                          % Maximum stiffness value (N/m)
st          = 0.001;                                                        % Simulation time step (s)
tDesMax     = 0;                                                            % Initialize maximum time duration

% Choose optimization algorithm (uncomment desired algorithm)
    algorithm = NLOPT_GN_CRS2_LM;                                       % Controlled Random Search (CRS) with local mutation
%     algorithm = NLOPT_GN_DIRECT;                                        % DIRECT (DIviding RECTangles)
%     algorithm = NLOPT_GN_DIRECT_NOSCAL;                                 % DIRECT, unscaled variant
%     algorithm = NLOPT_GN_DIRECT_L;                                      % DIRECT-L (locally biased variant of DIRECT)
%     algorithm = NLOPT_GN_DIRECT_L_RAND;                                 % DIRECT-L, slightly randomized variant
%     algorithm = NLOPT_GN_ISRES;                                         % ISRES (Improved Stochastic Ranking Evolution Strategy)
%     algorithm = NLOPT_GN_ESCH;                                          % ESCH (evolutionary algorithm)
    
% Run main program for selected experimental block and simulation settings
% Get information about experimental block
[subjNum, subjStr, trialDate, trialStr, blockStr, ~, ~,invalidTrials] = ...
    blockDictionary(blockNum);
disp('Experimental Block: ')
disp([subjStr,trialDate,trialStr,blockStr])

disp(subjNum)
    
% SAVEPEAKS
% Create structure to hold the peak data for this block of trials
% peakArray = cell(50,1);

% Number of rows and columns for plotting multiple figures in 1 window.
% For running an entire experimental block, choose 5 rows and 10 columns.
plotRows    = 5;
plotCols    = 10;

% Determine if current block is block 3 or block 4
block3or4 = blockNumFunc(blockNum);

if fitMethod == "fitToAverage"
    plotInd = true;
end

% To plot an individual trial
if plotInd
   % Choose specific trial number to plot
   trialNumber = 25;     % Index trial number starting from 1
   numStart = trialNumber - 1;
   numEnd = numStart;
   if fitMethod == "fitToAverage"
      numStart = 40;
      numEnd = 49;
      trialNumber = 'NaN';
   end
end

numTrials   = numEnd - numStart + 1;

% Timing and amplitude error bounds for input shaping impulses:
if timeMod == true
    shift   = 0.200;                                                        % Input shaping execution error time bound (s)
else
    shift   = 0.000;
end
ampError    = 0.00;                                                         % Max amplitude error (as +/- percentage of nominal)

% Set optimization intial values. Useful to modify if you just want to do one
% simulation with specific pre-determined values. Otherwise, start from
% B=10, K=100, and tdelay = delayMax/2
for initChoose = 0
    
    if initChoose == 0
        % Initial values at midpoint of parameter ranges
        bInit = (Bmax+Bmin)/2;
        kInit = (Kmax+Kmin)/2;
        tdelayInit = 2;
        
    elseif initChoose == 1
        % Typical intial values
        bInit = 10;
        kInit = 100;
        tdelayInit = 0;
        tdelayInit = 2*tdelayInit/delayMax+1;                                         % Scale to range between 1 and 3

    elseif initChoose == 2
        % Very stiff intial values
        bInit = 0;
        kInit = 750;
        tdelayInit = 0.125;
        tdelayInit = 2*tdelayInit/delayMax+1;                                         % Scale to range between 1 and 3

    elseif initChoose == 3
        % Very compliant intial values
        bInit = 50;
        kInit = 1;
        tdelayInit = 0.125;
        tdelayInit = 2*tdelayInit/delayMax+1;                                         % Scale to range between 1 and 3

    elseif initChoose == 4
        % Custom values
        bInit = 23.74;
        kInit = 357.5;
        tdelayInit = 0.250;
        tdelayInit = 2*tdelayInit/delayMax+1;                                         % Scale to range between 1 and 3
    end
end
    
% Define variable weights
wPos    = weights(1);
wTheta  = weights(2);
wVel    = weights(3);
wOmega  = weights(4);
wAcc    = weights(5);
wAlpha  = weights(6);

% Initialize table that will hold trial number, objective function and
% optimization parameter values.
% Note that row size is initialized to be at least the same as number of
% paramTable rows so there is space to display all optimization parameters.
numParams   = 15;
numRows = max(numTrials,numParams);
if optimizationType == "input shaping 4 impulse"
    numCols = 18;
elseif optimizationType == "input shaping 2 impulse impedance"
    numCols = 14;
elseif optimizationType == "input shaping 2 impulse no impedance"
    numCols = 12;
elseif optimizationType == "submovement"
    numCols = 11;
end
% returnArray = zeros(numRows,numCols);
returnArray = NaN(numRows,numCols);
    
% Initialize figure windows
if plotInd
    figure('Name','Cart velocity','units','normalized','outerposition',[0 0 1 1]);
    cv = get(gcf,'Number');
else
    cp = figure('Name','Cart position','units','normalized','outerposition',[0 0 1 1]);
    set(cp,'DefaultFigureVisible','off');

    bp = figure('Name','Ball angle','units','normalized','outerposition',[0 0 1 1]);
    set(bp,'DefaultFigureVisible','off');

    cv = figure('Name','Cart velocity','units','normalized','outerposition',[0 0 1 1]);
    set(cv,'DefaultFigureVisible','off');

    bv = figure('Name','Ball angular velocity','units','normalized','outerposition',[0 0 1 1]);
    set(bv,'DefaultFigureVisible','off');

    ca = figure('Name','Cart acceleration','units','normalized','outerposition',[0 0 1 1]);
    set(ca,'DefaultFigureVisible','off');

    ba = figure('Name','Ball angular acceleration','units','normalized','outerposition',[0 0 1 1]);
    set(ba,'DefaultFigureVisible','off');
end

% Initialize persistent variable array holding evaluation count for all trials
% G = globalData('evalCount', zeros(1,50));
    
% Create folder to hold this trial's data and subfolder to hold best-fit
% profiles
nameStr = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'));
topFolderStr = strcat("data/", modelStr, "/");
folderStr = strcat("data/",modelStr,"/",subjNum,"_B",block3or4,"/");
mkdir(folderStr, "best-fit profiles");

% Load start and stop indices of trials for this block
load(strcat("trim times/",blockStr),'Expression1')

% Initialize arrays to hold best-fit results after optimization
xopt_array = cell(50,2);

% Initialize stopping tolerance
relative_stop = 0.00001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLOPT FUNCTION
function [xopt, fmin] = NLopt(tc, tdes, xEnd, vStart, vEnd, pendIndex,...
        pos, vel, acc, theta, omega, alpha, num, algorithm, setMaxEval)
    
    % DEBUG
%     disp("inside NLopt")

    % Set optimization algorithm
    opt.algorithm = algorithm;
    
    % Define lower and upper parameter bounds
    if optimizationType == "input shaping 4 impulse"
        
        % Parameter array: [b, k, p, q, r, s, tdelay, fa11, fa12, fa21, fa22]
        if ampMod == true
            opt.lower_bounds = [Bmin, Kmin, 1, 1, 1, 1, 1, 1-ampError,...
                                1-ampError, 1-ampError, 1-ampError];
            opt.upper_bounds = [Bmax, Kmax, 3, 3, 3, 3, 3, 1+ampError,...
                                1+ampError, 1+ampError, 1+ampError];
                            
           % Initial guess: b=10, k=100, no timing or amplitude execution error
            init_guess = [bInit, kInit, 2, 2, 2, 2, tdelayInit, 1, 1, 1, 1];  
        else
            opt.lower_bounds = [Bmin, Kmin, 1, 1, 1, 1, 1];
            opt.upper_bounds = [Bmax, Kmax, 3, 3, 3, 3, 3];
            
            % Initial guess: b=10, k=100, no timing or amplitude execution error
            init_guess = [bInit, kInit, 2, 2, 2, 2, tdelayInit];    
        end

    elseif optimizationType == "input shaping 2 impulse impedance"
        
        % Parameter array: [b, k, p, q, tdelay, fa1, fa2]
        if ampMod == true
            opt.lower_bounds = [Bmin, Kmin, 1, 1, 1, 1-ampError, 1-ampError];
            opt.upper_bounds = [Bmax, Kmax, 3, 3, 3, 1+ampError, 1+ampError];
            
            % Initial guess: b=10, k=100, no timing or amplitude execution error
            init_guess = [bInit, kInit, 2, 2, tdelayInit, 1, 1]; 
        else
            opt.lower_bounds = [Bmin, Kmin, 1, 1, 1];
            opt.upper_bounds = [Bmax, Kmax, 3, 3, 3];
            
            % Initial guess: b=10, k=100, no timing or amplitude execution error
            init_guess = [bInit, kInit, 2, 2, tdelayInit];
        end
        
    elseif optimizationType == "input shaping 2 impulse no impedance"
        
        % Parameter array: [p, q, tdelay, fa1, fa2]
        if ampMod == true
            opt.lower_bounds = [1, 1, 1, 1-ampError, 1-ampError];                                
            opt.upper_bounds = [3, 3, 3, 1+ampError, 1+ampError];
            
            % Initial guess: no timing execution error, no amplitude error
            init_guess = [2, 2, tdelayInit, 1, 1];
        else
            opt.lower_bounds = [1, 1, 1];
            opt.upper_bounds = [3, 3, 3];
            
            % Initial guess: no timing execution error, no amplitude error
            init_guess = [2, 2, tdelayInit];
        end
        
    elseif optimizationType == "submovement"
        % Parameter array: [b, k, D1, tf1, ti2, tdelay]
        minSubDist = 0.05;                                                  % Minimum submovement distance (m)
        minSubDur = 0.200;                                                  % Minimum submovement duration (s)
        opt.lower_bounds = [Bmin, Kmin,    minSubDist,  minSubDur,        0,       1];
        opt.upper_bounds = [Bmax, Kmax, xEnd-minSubDist,   tdes,   tdes-minSubDur, 3];
        
        % Initial guess: b=10, k=100, equal amplitudes, equal durations, no timing delay error
        % Parameter array: [b, k, D1, tf1, ti2, tdelay]
        init_guess = [bInit, kInit, xEnd/2, tdes/2, tdes/2, tdelayInit];
    end

    % Objective function to be minimized
    opt.min_objective = @(x) objFunc(x, optimizationType,ampMod,...
        forwardF,ver,impedance,tc,tdes,delayMin,delayMax,xEnd,vStart,vEnd,...
        st,shift,simVersion,pendIndex,objective,pos,vel,...
        acc,theta,omega,alpha,weights,blockStr,blockNum,subjNum,num,tDesMax,...
        fitMethod);

    % Stopping criteria: stop optimization when an evaluation step
    % changes every component of x by less than xtol_rel multiplied by the
    % absolute value of that component of x (e.g., xtol_rel = .01 means the
    % optimization will stop when a step changes all parameters by less
    % than 1%)
    opt.xtol_rel = relative_stop;

    % Maximum number of evaluations before stopping the algorithm
    opt.maxeval = setMaxEval;

    % Run NLopt algorithm
    % @xopt: optimal values of the optimization hyperparameters
    % @fmin: minimum objective function value
    [xopt, fmin, ~] = nlopt_optimize(opt, init_guess);
    
%     G = globalData();
%     evals = G.evalCount;

%     % DEBUG
%     disp(evals)
% 
%     evalNum = evals(num+1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function handle for parallel version:
NLoptHandle = @NLopt;

% If the method is fitToAverage, first takes the average of all
% experimental trials in a single block. Then uses this average experimental
% trial as the single trial that the optimization algorithm fits to.
if fitMethod == "fitToAverage"
    
%     stNorm = st;
%     
%     % Initialize variables for longest & shortest durations of group
%     longestDuration = 0;
%     shortestDuration = 100;
%     for num = numStart:numEnd
%         numStr = num2str(num);
%         
%         % Construct file name and load file of one experimental trial
%         fileStr = strcat(subjStr,trialDate,trialStr,numStr);
%         load(fileStr,'pos','theta','vel','omega','acc','alpha','t')
%         
%         % Load start and stop indices of corresponding trial. Note that
%         % trials names are indexed from 0, so add 1 to access correct row
%         % in array
%         start   = Expression1(num+1,1);
%         stop    = Expression1(num+1,2);
%         
%         % DEBUG
% %         disp('Trial #:')
% %         disp('Trial duration (pre-computation):')
%         duration = t(stop)-t(start);
%         if duration > longestDuration
%            longestDuration = duration;
%         end
%         
%         if duration < shortestDuration
%             shortestDuration = duration;
%         end
%     end
%     
%     longestDuration = round(longestDuration,-log10(st));
% %     shortestDuration
%     normLen = int16(longestDuration/stNorm);
%     
%     % Initialize arrays to hold all experimental data
%     posArray    = zeros(numTrials,normLen);
%     thetaArray  = zeros(numTrials,normLen);
%     velArray    = zeros(numTrials,normLen);
%     omegaArray  = zeros(numTrials,normLen);
%     accArray    = zeros(numTrials,normLen);
%     alphaArray  = zeros(numTrials,normLen);
%     tdesArray   = zeros(numTrials,1);
%     
%     % Take average of selected range of experimental trials
%     for num = numStart:numEnd
%         numStr = num2str(num);
%         arrayNum = num-numStart+1;
% 
%         % Skip invalid trials
%         if ismember(num,invalidTrials)
%             continue
%         end
%     
%         % Construct file name and load file of one experimental trial
%         fileStr = strcat(subjStr,trialDate,trialStr,numStr);
%         load(fileStr,'pos','theta','vel','omega','acc','alpha','t')
% 
%         % Load start index of corresponding trial. Note that trials are
%         % indexed from 0 so add 1 to access correct row in array
%         start   = Expression1(num+1,1);
% 
%         % Trim experimental data at point where duration first exceeds
%         % longest trial duration. Function also evenly spaces data
%         % according to chosen time step.
%         stop = find(t > longestDuration + t(start),1);
%         if isempty(stop)
%             stop = length(t);
%         end
%         [pos,theta,vel,omega,acc,alpha,tdes,tc] = trimData(pos,theta,vel,...
%             omega,acc,alpha,t,st,start,stop);
% 
%         % Convert angles back to radians and position back to -0.125
%         theta = theta*(2*pi)/360;
%         omega = omega*(2*pi)/360;
%         alpha = alpha*(2*pi)/360;
%         pos = pos - 0.125;
% 
%         % Trim data again, this time to the exact longest duration
%         start = 1;
%         stop = int16(longestDuration/st);
%         
%         [pos,theta,vel,omega,acc,alpha,tdes,tc] = trimData(pos,theta,vel,...
%             omega,acc,alpha,tc,st,start,stop);
%         
%         % Add normalized profiles to arrays with all trials from group
%         % being averaged
%         posArray(arrayNum,:)      = pos;
%         thetaArray(arrayNum,:)    = theta;
%         velArray(arrayNum,:)      = vel;
%         omegaArray(arrayNum,:)    = omega;
%         accArray(arrayNum,:)      = acc;
%         alphaArray(arrayNum,:)    = alpha;
%         tdesArray(arrayNum,:)     = tdes;
%         
%     end
%     
%     % Delete skipped trials from arrays. Filter them out by finding the
%     % trials whose final position is 0.
%     [zeroRow, ~] = find(~posArray(:,end));
%     posArray(zeroRow,:)     = [];
% %     thetaArray(zeroRow,:)   = [];
%     velArray(zeroRow,:)     = [];
% %     omegaArray(zeroRow,:)   = [];
% %     accArray(zeroRow,:)     = [];
% %     alphaArray(zeroRow,:)   = [];
%     
%     % Take average of all trials
%     posAvg      = mean(posArray);
% %     thetaAvg    = mean(thetaArray);
%     velAvg      = mean(velArray);
% %     omegaAvg    = mean(omegaArray);
% %     accAvg      = mean(accArray);
% %     alphaAvg    = mean(alphaArray);
% 
%     % Define cart end position, start and end velocities, and pendulum
%     % release time from experimental data
%     xEnd        = posAvg(end);
%     vStart      = velAvg(1);
%     vEnd        = velAvg(end);
%     pendIndex   = find(velAvg > 0.1, 1);                                    % The first time index that cart velocity surpasses 0.1 m/s
% 
%     % Extend time vector to maximum delayed size. All time vectors will
%     % be the same length, but shorter time durations will be padded.
%     % DEBUG
% %     disp('Minimum time duration:')
%     tdes = shortestDuration;
% %     disp('tdes max:')
%     tDesMax = tdes + delayMax;
%     tc = 0:stNorm:tDesMax;
% 
%     % Initialize variables shared by objFunc.m objective function
%     G = globalData('evalCount', 0);
%     evalCount = G.evalCount;
% 
%     %%%%% Run optimization. Fit values to selected experimental trial. %%%%
%     [algorithm, lb, ub, maxeval] = NLopt();
% 
%     % Print trial number and total evaluation count
%     disp(['Trial: ',num2str(num)])
%     disp(['Number of evaluations: ', num2str(evalCount)])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting here is the standard method of finding a best-fit simulated
% trial for each individual experimental trial.
else

    % Slice start and stop indices into separate columns; this to avoid
    % making Expression1 a broadcast variable in the parfor loop
    starts = Expression1(:,1);
    stops = Expression1(:,2);

    % Loop through trials, fitting values for each one
    % Parallel version:
    for num = numStart:numEnd
%     parfor num = numStart:numEnd
        numStr = num2str(num);
        % DEBUG
%         disp(["Parfor num: ", numStr])

        % Construct file name and load file of one experimental trial
        fileStr = strcat(subjStr,trialDate,trialStr,numStr);
        
        % Parallel version:
        structTrial = load(fileStr,'pos','theta','vel','omega','acc','alpha','t');
        pos = structTrial.pos;
        theta = structTrial.theta;
        vel = structTrial.vel;
        omega = structTrial.omega;
        acc = structTrial.acc;
        alpha = structTrial.alpha;
        t = structTrial.t;
        
        % Load start and stop indices of corresponding trial. Note that trials
        % are indexed from 0, so add 1 to access correct row in array.
        start = starts(num+1);
        stop = stops(num+1);
%         start   = Expression1(num+1,1);
%         stop    = Expression1(num+1,2);

        % Trim experimental data
        [pos,theta,vel,omega,acc,alpha,tdes,~] = ...
            trimData(pos,theta,vel,omega,acc,alpha,t,st,start,stop);

        % Define cart final position, initial and final velocities, and
        % pendulum release time from experimental data
        xEnd = pos(end);
        vStart = vel(1);
        vEnd = vel(end);
        pendIndex = find(vel > 0.1, 1);                                     % The first time index that cart velocity surpasses 0.1 m/s

        % Extend time vector to maximum delayed size. All time vectors will
        % be the same length, so shorter trials will be padded.
        tDesMax = tdes + delayMax;
        tc=0:st:tDesMax;

        %%% Run optimization. Fit values to selected experimental trial. %%%
        [xopt, fmin] = feval(NLoptHandle,tc,tdes,xEnd,vStart,vEnd,...
            pendIndex,pos,vel,acc,theta,omega,alpha,num,algorithm,setMaxEval);

        % Store best-fit parameters and objective function value in array 
        % to be accessed outside parfor-loop
        xopt_array(num+1,:) = {xopt, fmin};

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output results of optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for num = numStart:numEnd

    % DEBUG
%     disp("Outputting results")

    xopt = xopt_array{num+1,1};
    fmin = xopt_array{num+1,2};

    % Extract best-fit parameters from optimization output
    optB = xopt(1);
    optK = xopt(2);
    if optimizationType == "input shaping 4 impulse"
        optP        = xopt(3);
        optQ        = xopt(4);
        optR        = xopt(5);
        optS        = xopt(6);
        optTcount   = xopt(7);
        
        if ampMod == true
            optA11  = xopt(8);
            optA12  = xopt(9);
            optA21  = xopt(10);
            optA22  = xopt(11);
        else
            optA11  = 1;
            optA12  = 1;
            optA21  = 1;
            optA22  = 1;
        end
        
    elseif optimizationType == "input shaping 2 impulse impedance"
        optP        = xopt(3);
        optQ        = xopt(4);
        optTcount   = xopt(5);
        
        if ampMod == true
            optA1   = xopt(6);
            optA2   = xopt(7);
        else
            optA1   = 1;
            optA2   = 1;
        end
    elseif optimizationType == "input shaping 2 impulse no impedance"
        optP        = xopt(1);                                              
        optQ        = xopt(2);
        optTcount   = xopt(3);
        
        if ampMod == true
            optA1   = xopt(4);
            optA2   = xopt(5);
        else
            optA1   = 1;
            optA2   = 1;
        end
    elseif optimizationType == "submovement"
        optD        = xopt(3);
        optt1       = xopt(4);
        optt2       = xopt(5);
        optTcount   = xopt(6);
    end
    %%%%%%

    % Construct file name and load file of one experimental trial
    numStr = num2str(num);
    fileStr = strcat(subjStr, trialDate, trialStr, numStr);
    
    load(fileStr, 'pos', 'theta', 'vel', 'omega', 'acc', 'alpha', 't');
%     structTrial = load(fileStr,'pos','theta','vel','omega','acc','alpha','t');
%     pos = structTrial.pos;
%     theta = structTrial.theta;
%     vel = structTrial.vel;
%     omega = structTrial.omega;
%     acc = structTrial.acc;
%     alpha = structTrial.alpha;
%     t = structTrial.t;

    % Load start and stop indices of corresponding trial. Note that trials
    % are indexed from 0, so add 1 to access correct row in array.
    start = Expression1(num+1, 1);
    stop = Expression1(num+1, 2);

    % Trim experimental data. trimData takes angles in radians and returns
    % them in degrees.
    [pos, theta, vel, omega, acc, alpha, tdes, ~] = ...
        trimData(pos, theta, vel, omega, acc, alpha, t, st, start, stop);

    % Define cart final position, initial and final velocities, and
    % pendulum release time from experimental data
    xEnd = pos(end);
    vStart = vel(1);
    vEnd = vel(end);
    pendIndex = find(vel > 0.1, 1); % first time index when cart velocity exceeds 0.1 m/s

    % Extend time vector to maximum delayed size. All time vectors will
    % be the same length, so shorter trials will be padded.
    tDesMax = tdes + delayMax;
    tc = 0:st:tDesMax;

    % Translate optimization time parameters to actual times in seconds
    optTdelay = 0.5 * (optTcount - 1) * delayMax; % Range from delayMin to delayMax
%     optTdelay   = (optTcount-2)*delayMax;     % Range from -delayMax to +delayMax (MARK FOR DELETION)
    tdessim = tdes + optTdelay;
    if optimizationType == "input shaping 4 impulse"
        optShift = shift * [optP - 2, optQ - 2, optR - 2, optS - 2];
    elseif optimizationType == "input shaping 2 impulse impedance" ||...
            optimizationType == "input shaping 2 impulse no impedance"
        optShift = shift * [optP - 2, optQ - 2];
    end
    
    % Create system using optimal hyperparameters
    [sys, sysRigid, Td, Td1, Td2, zeta, zeta1, zeta2, overdamped] = ...
        sysCreate(optB, optK, forwardF, "full", impedance, printSys);

    if printSys
        % Output damping ratios and damped natural periods
        optB
        optK
        zeta1
        zeta2
        Td1
        Td2
    end    

    % Simulate system with optimized parameters
    if optimizationType == "input shaping 4 impulse"
        modes = 2;
        [output, vc, sub1End, sub2Start, dur_corr] = simInputShape(optB, ...
            optK, ver, sys, sysRigid, Td1, Td2, zeta1, zeta2, tdes, ...
            tdessim, xEnd, vStart, vEnd, optA11, optA12, optA21, optA22, ...
            optP, optQ, optR, optS, st, shift, forwardF, simVersion, ...
            modes, pendIndex, fitMethod);
        
    elseif optimizationType == "input shaping 2 impulse impedance" || ...
            optimizationType == "input shaping 2 impulse no impedance"
        
        modes   = 1;
        Td2     = 0;
        zeta2   = 0;
        fa21    = 0;
        fa22    = 0;
        r       = 0;
        s       = 0;
        [output,vc,sub1End,sub2Start,dur_corr] = simInputShape(optB,optK,ver,sys,...
            sysRigid,Td,Td2,zeta,zeta2,tdes,tdessim,xEnd,vStart,vEnd,...
            optA1,optA2,fa21,fa22,optP,optQ,r,s,st,shift,forwardF,...
            simVersion,modes,pendIndex,fitMethod);

    elseif optimizationType == "submovement"
        [output, vc] = simSubmovements(sys,sysRigid,tdes,tdessim,st,xEnd,vStart,...
            optD,optt1,optt2,pendIndex,simVersion,forwardF,impedance);
    end
    
    % If using feature-based objective function, find location of simulated
    % peaks
    if objective == "features"
       objCalc(objective,output,pos,theta,vel,omega,acc,alpha,weights,...
           blockStr,blockNum,subjNum,num);
       G = globalData();
       velChange_exp = G.velChange_exp;
       velPeaks_exp = G.velPeaks_exp;
       velChange_sim = G.velChange_sim;
       velPeaks_sim = G.velPeaks_sim;

       % SAVEPEAKS
       % Add peaks to cell array
%        peakArray{num+1,1} = [velChange_exp; velPeaks_exp];

    % If just plotting peaks, without using for objective function
    elseif plotPeaks
        vel_sim = output(:,3);
        acc_sim = output(:,5); 
        [velChange_sim, velPeaks_sim] = peakCalc(vel_sim, acc_sim);
    end
    
    % Rename outputs for ease of use
    pos_sim     = output(:, 1);
    theta_sim   = rad2deg(output(:, 2));
    vel_sim     = output(:, 3);
    omega_sim   = rad2deg(output(:, 4));
    acc_sim     = output(:, 5);
    alpha_sim   = rad2deg(output(:, 6));
    
    % Handle any size difference due to rounding error
    lensim = length(vel_sim);
    lenexp = length(pos);
    lentc  = length(tc);
    if lensim ~= lenexp
        smaller = min(lensim,lenexp);

        % To calculate RMSE only during times when both profiles are
        % defined, shorten the longer array to the shorter one's length
        if lensim > lenexp
            vel_sim_RMSE    = vel_sim(1:smaller);
            vel_RMSE        = vel;
        else
            vel_RMSE        = vel(1:smaller);
            vel_sim_RMSE    = vel_sim;
        end
        
    else
        vel_RMSE            = vel;
        vel_sim_RMSE        = vel_sim;
    end
    
    % Make time vector same length as simulated data
    if lensim ~= lentc
        tc = 0:st:(lensim-1)*st;
    end
    
    % Calculate cart velocity RMSE
    dvel        = vel_sim_RMSE - vel_RMSE';
    dvel        = dvel/max(vel);
    squareVel   = dvel.^2;
    vrmse       = sqrt((1/length(squareVel))*sum(squareVel));
    
    % Insert trial number, fmin and optimized parameter values into array
    returnRow = num-numStart+1;
    trial = num+1;
    if fitMethod == "fitToAverage"
        returnRow = 1;
        trial = NaN; 
    end

    % Fill array with best-fit optimization parameters
    if optimizationType == "input shaping 4 impulse"
        returnArray(returnRow,:) = [trial, fmin, vrmse, optB, optK,...
            optShift,optTdelay, optA11, optA12, optA21, optA11, sub1End,...
            sub2Start, 0, relative_stop];
    elseif optimizationType == "input shaping 2 impulse impedance"
        returnArray(returnRow,:) = [trial, fmin, vrmse, optB, optK,...
            optShift, optTdelay, optA1, optA2, sub1End, sub2Start,...
            0, relative_stop];
    elseif optimizationType == "input shaping 2 impulse no impedance"
        returnArray(returnRow,:) = [trial, fmin, vrmse, optShift,...
            optTdelay, optA1, optA2, sub1End, sub2Start, 0,...
            relative_stop];
    elseif optimizationType == "submovement"
        returnArray(returnRow,:) = [trial, fmin, vrmse, optB, optK, optD,...
            optt1, optt2, optTdelay, 0, relative_stop];
    end
    
    % Save best-fit output to a spreadsheet    
    dur_corr_vector = zeros(length(vel_sim),1);
    dur_corr_vector(1) = dur_corr;
    outputArray = [tc', pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim, dur_corr_vector];
    profileVariables = {'Time [s]', 'Cart Pos [m]', 'Ball Angle [deg]',...
        'Cart Vel [m/s]','Ball Ang Vel [deg/s]', 'Cart Acc [m/s^2]',...
        'Ball Ang Acc [deg/s^2]', 'Pre-Trim Duration [s]'};
    ProfileT = array2table(outputArray,'VariableNames',profileVariables);
    trialNumStr = num2str(returnRow);
    fileName = strcat(folderStr, "/best-fit profiles/", trialNumStr, ".xlsx");
    writetable(ProfileT,fileName);
    
    %%%% Check if "best" solution is actually one of filtered cases
    % If cart velocity solution is ~0 velocity
    if mean(vel_sim)<0.1      
        disp('Chosen trajectory is 0 velocity');

    % If cart velocity exceeds saturation limit
    elseif max(vel_sim) > 0.5
        disp('Chosen trajectory over max velocity');

    % If cart velocity exceeds saturation minimum
    elseif min(vel_sim) < -0.075
        disp('Chosen trajectory under min velocity');

    % If ball angle exceeds saturation limit
    elseif max(theta_sim) > 50
        disp('Chosen trajectory over max ball angle');

    % If ball angular velocity exceeds saturation limit
    elseif max(omega_sim) > 75
        disp('chosen trajectory over max ball velocity');

    % If cart acceleration exceeds saturation limit
    elseif max(acc_sim) > 5
        disp('Chosen trajectory over max cart acceleration');
        %{
    % If ball acceleration exceeds saturation limit
    elseif max(alpha_sim) > 15
        disp('Chosen trajectory over max ball acceleration')
        %}
    end
    
    % Define time vector for experimental profile
    te = 0:st:(lenexp-1)*st;
    
    % Subplot current trial in large figure window with multiple trials
    if plotInd
        % Plot Individual Velocity Trajectory
        labelFontSize = 24;
        titleFontSize = 28;
        axisFontSize = 16;
        lineWidth = 3;
        colGreen = [0, 153, 51]/255;
        colOrange = [255, 153, 0]/255;
        colBlue = [0, 102, 255]/255;
        colGrey = [153, 153, 153]/255;
        colDarkGrey = [51, 51, 51]/255;
        defaultBlue = [0, 0.4470, 0.7410];
        defaultOrange = [0.8500, 0.3250, 0.0980];

        % For characteristic trials:
        % Cart velocity (vel): plot experimental versus simulation
        set(groot,'CurrentFigure',cv);
        plot(te,vel,'LineWidth',lineWidth,'Color',defaultBlue)
        hold on;
        if plotSim      % Plot simulated profile
            plot(tc,vel_sim,'LineWidth',lineWidth,'Color',defaultOrange)
        end
        if plotPeaks   % Plot local velocity max & min
           plot(velChange_exp*st,velPeaks_exp,'Marker','*','MarkerSize',12,...
               'Color','magenta','LineWidth',2,'LineStyle','none')
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Time (s)','FontSize',labelFontSize)
        ylabel('Cart Velocity (m/s)','FontSize',labelFontSize)
        subjNum = erase(subjNum,'S');
        if fitMethod ~= "fitToAverage"
            title(['Subject ', subjNum, ', Block ', block3or4 ', Trial ',...
                num2str(trialNumber)],'FontSize',titleFontSize);
        end
        
        % For flow diagram: simulated cart velocity only
        %{
        labelFontSize = 18;
        titleFontSize = 20;
        axisFontSize = 16;
        lineWidth = 2;
        
        figure();
        plot(tc,vel_sim,'LineWidth',lineWidth,'Color',defaultOrange)
        set(gca,'FontSize',axisFontSize)
        xlabel('Time (s)','FontSize',labelFontSize)
        ylabel('Cart Velocity (m/s)','FontSize',labelFontSize)
        title('Output Cart Velocity', 'FontSize',titleFontSize)
        
        % Trim simulated profile
        if length(vel_sim) > length(vel)
            vel_sim     = vel_sim(1:length(vel));
            tc          = tc(1:length(vel));
        end
        
        % For flow diagram: both simulated & experimental cart velocity
        figure();
        plot(te,vel,'LineWidth',lineWidth,'Color',defaultBlue)
        hold on;
        plot(tc,vel_sim,'LineWidth',lineWidth,'Color',defaultOrange)
        set(gca,'FontSize',axisFontSize)
        xlabel('Time (s)','FontSize',labelFontSize)
        ylabel('Cart Velocity (m/s)','FontSize',labelFontSize)
        title('Experimental and Simulated Cart Velocities','FontSize',titleFontSize)
        %}
        
    else
        
        % Before sending data to plot, make sure all arrays are row vectors
        sz_te = size(te);
        sz_tc = size(tc);
        sz_pos = size(pos);
        sz_pos_sim = size(pos_sim);

        if sz_te(1) ~= 1
            te = te';
        end

        if sz_tc(1) ~= 1
            tc = tc';
        end

        if sz_pos(1) ~=1
            pos = pos';
            % ADD OTHER VARIABLES HERE ONCE IT WORKS FOR POS
        end
        
        if sz_pos_sim(1) ~= 1
            pos_sim = pos_sim';
        end

        % Trial number as string for plot titles
        trialCaption = num2str(num+1);

        % Cart position (pos): plot experimental versus simulation
%         posData = cat(2, plotRows, plotCols, num, numStart, length(te), length(tc), te, pos, tc, pos_sim);
%         send(D, posData)

        set(groot,'CurrentFigure',cp);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,pos,'LineWidth',2)
        hold on;
        if plotSim
            plot(tc,pos_sim,'LineWidth',2)
        end
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('cart position (m)')

        % Ball angle (theta): plot experimental versus simulation
        set(groot,'CurrentFigure',bp);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,theta,'LineWidth',2)
        hold on;
        if plotSim
            plot(tc,theta_sim,'LineWidth',2)
        end
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('ball angle (deg)')

        % Cart velocity (vel): plot experimental versus simulation
        set(groot,'CurrentFigure',cv);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,vel,'LineWidth',2)
        hold on;
        if plotSim && ~printDes
            plot(tc,vel_sim,'LineWidth',2)
%             if plotPeaks && (objective == "features")   % Plot simulated local velocity max & min
            if plotPeaks   % Plot simulated local velocity max & min
                plot(velChange_sim*st,velPeaks_sim,'Marker','s','MarkerSize',12,...
               'Color','cyan','LineWidth',1.5,'LineStyle','none')
            end
        end
        if plotPeaks && (objective == "features")   % Plot experimental local velocity max & min
            plot(velChange_exp*st,velPeaks_exp,'Marker','s','MarkerSize',12,...
               'Color','magenta','LineWidth',1.5,'LineStyle','none')
        end
        if printDes
            % Plot individual input velocity submovements
            [vcRows, vcCols] = size(vc);
            td = 0:st:(vcCols-1)*st;
            for row = 1:vcRows
                plot(td,vc(row,:),'LineWidth',1.25,'Color',[0.9290, 0.6940, 0.1250])   
            end
            % Plot simulated output velocity
            plot(tc, vel_sim, 'LineWidth', 2,'LineStyle','--','Color',[0.8500, 0.3250, 0.0980]);
        end
        title(['Trial' ' ' trialCaption]);
        xlabel('Time (s)')
        ylabel('Cart Velocity (m/s)')

        % Ball angular velocity (omega): plot experimental versus simulation
        set(groot,'CurrentFigure',bv);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,omega,'LineWidth',2)
        hold on;
        if plotSim
            plot(tc,omega_sim,'LineWidth',2)
        end
        trial = num2str(num+1);
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('ball velocity (deg/s)')

        % Cart acceleration (acc): plot experimental versus simulation
        set(groot,'CurrentFigure',ca);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,acc,'LineWidth',2)
        hold on;
        if plotSim
            plot(tc,acc_sim,'LineWidth',2)
        end
        trial = num2str(num+1);
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('cart acceleration (m/s^2)')

        % Ball acceleration (alpha): plot experimental versus simulation
        set(groot,'CurrentFigure',ba);
        subplot(plotRows,plotCols,num-numStart+1)                               
        plot(te,alpha,'LineWidth',2)
        hold on;
        if plotSim
            plot(tc,alpha_sim,'LineWidth',2)
        end
        trial = num2str(num+1);
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('ball acceleration (deg/s^2)')
        
    end    

%     % Get number of evaluations after optimization
%     G = globalData();
%     evals = G.evalCount;
%     disp(evals)
% 
%     % Print trial number and total evaluation count
%     disp(['Trial: ',num2str(num+1)])
%     disp(['Number of evaluations: ', num2str(evals(num+1))])

end



% For determining what algorithm name the NLopt number corresponds to.
% disp('algorithm number: ')
% disp(num2str(algorithm))

% % Retrieve optimization info from persistent variable
% algorithm = G.algorithm;
% maxeval = G.maxeval;

% Look up NL opt algorithm name from double output
if algorithm == 0
    algorithm = "NLOPT_GN_DIRECT";
elseif algorithm == 1
    algorithm = "NLOPT_GN_DIRECT_L";
elseif algorithm == 2
    algorithm = "NLOPT_GN_DIRECT_L_RAND";
elseif algorithm == 3
    algorithm = "NLOPT_GN_DIRECT_NOSCAL";
elseif algorithm == 19
    algorithm = "NLOPT_GN_CRS2_LM";
elseif algorithm == 35
    algorithm = "NLOPT_GN_ISRES";
elseif algorithm == 42
    algorithm = "NLOPT_GN_ESCH";
else
    algorithm = "UNKNOWN";
end

% Convert optimization parameters to strings to include in output table
blocknameStr    = strcat("Block: ",subjStr,trialDate,trialStr);            
blocknameStr    = erase(blocknameStr,"_trial_"); 

algorithm       = strcat("Algorithm: ", algorithm);
Bbounds         = strcat("B bounds: ", convertCharsToStrings(num2str(Bmin)),...
                    " to ", convertCharsToStrings(num2str(Bmax)));
Kbounds         = strcat("K bounds: ", convertCharsToStrings(num2str(Kmin)),...
                    " to ", convertCharsToStrings(num2str(Kmax)));
durationtime    = strcat("Max duration error: ",...
                    convertCharsToStrings(num2str(delayMax*1000)), " ms");
setMaxEval         = strcat("Eval max: ",...
                    convertCharsToStrings(num2str(setMaxEval)));
wPosStr         = strcat("Cart position weight: ",...
                    convertCharsToStrings(num2str(wPos)));
wThetaStr       = strcat("Ball position weight: ",...
                    convertCharsToStrings(num2str(wTheta)));
wVelStr         = strcat("Cart velocity weight: ",...
                    convertCharsToStrings(num2str(wVel)));
wOmegaStr       = strcat("Ball velocity weight: ",...
                    convertCharsToStrings(num2str(wOmega)));
wAccStr         = strcat("Cart acceleration weight: ",...
                    convertCharsToStrings(num2str(wAcc)));
wAlphaStr       = strcat("Ball acceleration weight: ",...
                    convertCharsToStrings(num2str(wAlpha)));

if optimizationType == "input shaping 4 impulse" ||...
        optimizationType == "input shaping 2 impulse impedance" ||...
        optimizationType == "input shaping 2 impulse no impedance"
    impulsetime = strcat("Max impulse error: ", convertCharsToStrings(num2str(shift*1000)), " ms");
    amplitude   = strcat("Max amplitude error: ", convertCharsToStrings(num2str(100*ampError)), "%");
end

% Add current block name and optimization settings to table
paramTable          = table('Size',[numTrials,1],'VariableTypes',{'string'},... % Initialize 1-column table
                    'VariableNames',{'Parameters'});
paramTable{1,1}     = blocknameStr;
paramTable{3,1}     = algorithm;
paramTable{4,1}     = Bbounds;
paramTable{5,1}     = Kbounds;
paramTable{8,1}     = durationtime;
paramTable{9,1}     = setMaxEval;
paramTable{10,1}    = wPosStr;
paramTable{11,1}    = wThetaStr;
paramTable{12,1}    = wVelStr;
paramTable{13,1}    = wOmegaStr;
paramTable{14,1}    = wAccStr;
paramTable{15,1}    = wAlphaStr;

if optimizationType == "input shaping 4 impulse"
    paramTable{2,1} = "4-Impulse Input Shaping Fitting";
    paramTable{6,1} = impulsetime;
    paramTable{7,1} = amplitude;
elseif optimizationType == "input shaping 2 impulse impedance"
    paramTable{2,1} = "2-Impulse Input Shaping w/ Impedance Fitting";
    paramTable{6,1} = impulsetime;
    paramTable{7,1} = amplitude;
elseif optimizationType == "input shaping 2 impulse no impedance"
    paramTable{2,1} = "2-Impulse Input Shaping w/out Impedance Fitting";
    paramTable{4,1} = " ";                                                      % Overwrite B and K bounds with empty string
    paramTable{5,1} = " ";
    paramTable{6,1} = impulsetime;
    paramTable{7,1} = amplitude;
elseif optimizationType == "submovement"
    paramTable{2,1} = "Submovement Fitting";
    paramTable{6,1} = "";
end 

%%%% Save return array to a file
% Convert array to table
if optimizationType == "input shaping 4 impulse"
    variableNames = {'Trial','Fmin','VRMSE','B','K','Impulse1','Impulse2',...
                    'Impulse3','Impulse4','Tdelay','Amp11','Amp12','Amp21',...
                    'Amp22','Sub1End','Sub2Start','Evaluations','Stop Tolerance %'};
elseif optimizationType == "input shaping 2 impulse impedance"
    variableNames = {'Trial','Fmin','VRMSE','B','K','Impulse1','Impulse2',...
                    'Tdelay','Amp1','Amp2','Sub1End','Sub2Start',...
                    'Evaluations','Stop Tolerance %'};
elseif optimizationType == "input shaping 2 impulse no impedance"
    variableNames = {'Trial','Fmin','VRMSE','Impulse1','Impulse2',...
                    'Tdelay','Amp1','Amp2','Sub1End','Sub2Start',...
                    'Evaluations','Stop Tolerance %'};
elseif optimizationType == "submovement"
    variableNames = {'Trial','Fmin','VRMSE','B','K','Dist1',...
                    'Submovement 1 end','Submovement 2 start','Tdelay',...
                    'Evaluations','Stop Tolerance %'};
end

returnTable = array2table(returnArray,'VariableNames',variableNames);

% Concatenate tables for output
T = [paramTable, returnTable];

% Analyze data in array of optimized input parameters
if optimizationType == "input shaping 4 impulse"
    analysisColumns = 16;
    analysisTypes = {'string','double','double','double','double',...
                    'double','double','double','double','double',...
                    'double','double','double','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','B','K','Impulse1',...
                    'Impulse2','Impulse3','Impulse4','Amp11','Amp12',...
                    'Amp21','Amp22','Sub1End','Sub2Start','Tdelay'};
elseif optimizationType == "input shaping 2 impulse impedance"
    analysisColumns = 12;
    analysisTypes = {'string','double','double','double','double','double',...
                    'double','double','double','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','B','K','Impulse1',...
                    'Impulse2','Tdelay','Amp1','Amp2','Sub1End','Sub2Start'};
elseif optimizationType == "input shaping 2 impulse no impedance"
    analysisColumns = 10;
    analysisTypes = {'string','double','double','double','double','double'...
                    'double','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','Impulse1','Impulse2',...
                    'Tdelay','Amp1','Amp2','Sub1End','Sub2Start'};
elseif optimizationType == "submovement"
    analysisColumns = 9;
    analysisTypes = {'string','double','double','double','double',...
                    'double','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','B','K','Dist1','Sub1End',...
                    'Sub2Start','Tdelay'};
end

AnalysisT = table('Size',[size(returnTable,1), analysisColumns],...
            'VariableTypes',analysisTypes,'VariableNames',analysisNames);

% Populate analysis table with minimum, maximum, mean, median, and standard
% deviation of optimization parameters
AnalysisT{1,1} = "Min";
AnalysisT{2,1} = "Max";
AnalysisT{3,1} = "Mean";
AnalysisT{4,1} = "Median";
AnalysisT{5,1} = "Std Dev";

for i = 2:analysisColumns
    
    AnalysisT{1,i} = min(returnArray(:,i));
    AnalysisT{2,i} = max(returnArray(:,i));
    AnalysisT{3,i} = mean(returnArray(:,i));
    AnalysisT{4,i} = median(returnArray(:,i));
    AnalysisT{5,i} = std(returnArray(:,i));
    
end

% Save tables to excel file
% Note: nameStr and folderStr were already defined before trials were
% looped through 
fileName    = strcat(folderStr, nameStr, ".xlsx");
writetable(T,fileName);                                                     % Save paramTable and returnTable as excel file in new folder
writetable(AnalysisT,fileName,'Range','U1:BB6');                            % Add AnalysisT to excel file

% Save just data without headings to a .mat file
save(strcat(topFolderStr, num2str(blockNum), ".mat"), "returnArray");

% Convert strings to character arrays for plotting
subjStr     = strcat(char(subjStr),'\');
trialDate   = char(trialDate);

% Give title to and save figures
if plotInd
    set(groot,'CurrentFigure',cv);
    if fitMethod == "fitToAverage"
        figTitleStr = 'Average of Last 10 Hand Velocity Profiles and Best-Fit Simulation: Subject ';
        figTitle = [figTitleStr, subjNum, ' Block ', block3or4 ];
        sgtitle(figTitle,'FontSize',16)
    end
    saveas(cv, strcat(folderStr, "cart_vel_", nameStr, ".png"));
    saveas(cv, strcat(topFolderStr, "cart_vel_", nameStr, ".png"));         % Save velocity profile into top-level folder as well
else
%     set(groot,'CurrentFigure',cp);
%     set(cp,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Cart Position Profiles,' ' ' subjStr trialDate];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(cp,figTitle,'FontSize',16)
    
%     set(groot,'CurrentFigure',bp);
%     set(bp,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Ball Angle Profiles,' ' ' subjStr trialDate];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(bp,figTitle,'FontSize',16)

%     set(groot,'CurrentFigure',cv);
%     set(cv,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Cart Velocity Profiles: Subject ', subjNum, ' Block ', block3or4 ];
%     figTitle = ['Hand Velocity Minima & Maxima: Subject ', subjNum, ' Block ', block3or4 ];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(cv,figTitle,'FontSize',16)

%     set(groot,'CurrentFigure',bv);
%     set(bv,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Ball Angular Velocity Profiles,' ' ' subjStr trialDate];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(bv,figTitle,'FontSize',16)

%     set(groot,'CurrentFigure',ca);
%     set(ca,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Cart Acceleration Profiles,' ' ' subjStr trialDate];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(ca,figTitle,'FontSize',16)

%     set(groot,'CurrentFigure',ba);
%     set(ba,'DefaultFigureVisible','on');
    figTitle = ['Experimental vs Best Fit Ball Angular Acceleration Profiles,' ' ' subjStr trialDate];
%     sgtitle(figTitle,'FontSize',16)
    sgtitle(ba,figTitle,'FontSize',16)

    saveas(cp, strcat(folderStr, "cart_pos_", nameStr, ".png"))      
    saveas(bp, strcat(folderStr, "ball_angle_", nameStr, ".png"))      
    saveas(cv, strcat(folderStr, "cart_vel_", nameStr, ".png"))
    saveas(bv, strcat(folderStr, "ball_angular_vel_", nameStr, ".png"))
    saveas(ca, strcat(folderStr, "cart_acc_", nameStr, ".png"))
    saveas(ba, strcat(folderStr, "ball_angular_acc_", nameStr, ".png"))
end

% SAVEPEAKS
% Save cell array to a .mat file
% peakFileName = strcat("peak times/", blockStr);
% save(peakFileName, 'peakArray')

% Give focus back to command window when saving is complete
commandwindow

end

