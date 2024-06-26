function optimization(blockNum,optimizationType,modelStr,forwardF,intModel,...
        impedance,simVersion,weights,delayMin,delayMax,plotInd,printDes,...
        printSys,plotSim,plotPeaks,setMaxEval,numStart,numEnd)
% OPTIMIZATION
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
%
% * Physical plant
%   ** Including vs. not including hand impedance
%   ** Linear approximation vs. exact nonlinear model
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
% function optimizationCall.m automatically does. (Adds folders "data" and 
% "trim times" to the path, as well as "experimental data" and its
% subfolders.)

% Set optimization parameter values. B is hand damping in N*s/m, K is hand
% stiffness in N/m, st is simulation time step in seconds, and tDesMax is 
% initialized to 0. Note that constants for gravity, and physical system
% mass and length are defined in globalData.m)
Bmin = 0;
Bmax = 100;
Kmin = 0;
Kmax = 1000;
st = 0.001;
tDesMax = 0;

% Choose optimization algorithm (uncomment desired algorithm)
% Controlled Random Search (CRS) with local mutation:
algorithm = NLOPT_GN_CRS2_LM;
% DIRECT (DIviding RECTangles):
% algorithm = NLOPT_GN_DIRECT;
% DIRECT, unscaled variant:
% algorithm = NLOPT_GN_DIRECT_NOSCAL;
% DIRECT-L (locally biased variant of DIRECT):
% algorithm = NLOPT_GN_DIRECT_L;
% DIRECT-L, slightly randomized variant:
% algorithm = NLOPT_GN_DIRECT_L_RAND;
% ISRES (Improved Stochastic Ranking Evolution Strategy):
% algorithm = NLOPT_GN_ISRES;
% ESCH (evolutionary algorithm):
% algorithm = NLOPT_GN_ESCH;
    
% Run main program for selected experimental block and simulation settings
% Get information about experimental block
[subjNum, subjStr, trialDate, trialStr, blockStr,~,~,~] = ...
    blockDictionary(blockNum);
disp('Experimental Block: ')
disp([subjStr,trialDate,trialStr,blockStr])
disp(subjNum)

% Number of rows and columns for plotting multiple figures in 1 window --
% 5 rows by 10 columns to visualize the 50 trials in a single block
plotRows = 5;
plotCols = 10;

% Determine if current block is block 3 or block 4
block3or4 = blockNumFunc(blockNum);

% To plot an individual trial
if plotInd
   trialNumber = 25;     % Index trial number starting from 1
   numStart = trialNumber - 1;
   numEnd = numStart;
end

numTrials = numEnd - numStart + 1;

% Set optimization intial values. Useful to modify if running simulation
% with specific pre-determined values. Otherwise, start from B=10, K=100,
% and tdelay = delayMax/2
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
        tdelayInit = 2*tdelayInit/delayMax+1;

    elseif initChoose == 2
        % Very stiff intial values
        bInit = 0;
        kInit = 750;
        tdelayInit = 0.125;
        tdelayInit = 2*tdelayInit/delayMax+1;

    elseif initChoose == 3
        % Very compliant intial values
        bInit = 50;
        kInit = 1;
        tdelayInit = 0.125;
        tdelayInit = 2*tdelayInit/delayMax+1;

    elseif initChoose == 4
        % Custom values
        bInit = 23.74;
        kInit = 357.5;
        tdelayInit = 0.250;
        tdelayInit = 2*tdelayInit/delayMax+1;
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
% optimization parameter values
numParams   = 15;
numRows = max(numTrials,numParams);
if optimizationType == "input shaping 4 impulse" || ...
        optimizationType == "input shaping 2 impulse impedance"
    numCols = 6;
elseif optimizationType == "input shaping 2 impulse no impedance"
    numCols = 4;
end
returnArray = NaN(numRows,numCols);
    
% Initialize figure windows
if plotInd
    figure('Name','Cart velocity','units','normalized',...
        'outerposition',[0 0 1 1]);
    cv = get(gcf,'Number');
else
    cp = figure('Name','Cart position','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(cp,'DefaultFigureVisible','off');

    bp = figure('Name','Ball angle','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(bp,'DefaultFigureVisible','off');

    cv = figure('Name','Cart velocity','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(cv,'DefaultFigureVisible','off');

    bv = figure('Name','Ball angular velocity','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(bv,'DefaultFigureVisible','off');

    ca = figure('Name','Cart acceleration','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(ca,'DefaultFigureVisible','off');

    ba = figure('Name','Ball angular acceleration','units','normalized',...
        'outerposition',[0 0 1 1]);
    set(ba,'DefaultFigureVisible','off');
end
    
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
% NLopt Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xopt, fmin] = NLopt(tc, tdes, xEnd, vStart, pendIndex, pos, vel,...
        acc, theta, omega, alpha, algorithm, setMaxEval)

    % Set optimization algorithm
    opt.algorithm = algorithm;
    
    % Define parameter bounds and intial guess
    if optimizationType == "input shaping 4 impulse" || ...
            optimizationType == "input shaping 2 impulse impedance"
        % Parameter array: [b, k, tdelay]
        opt.lower_bounds = [Bmin, Kmin, 1];
        opt.upper_bounds = [Bmax, Kmax, 3];
        init_guess = [bInit, kInit, tdelayInit];
    elseif optimizationType == "input shaping 2 impulse no impedance"
        % Parameter array: [tdelay]
        opt.lower_bounds = [1];
        opt.upper_bounds = [3];
        init_guess = [tdelayInit];
    end

    % Objective function to be minimized
    opt.min_objective = @(x) objFunc(x, optimizationType,forwardF,intModel,...
        impedance,tc,tdes,delayMin,delayMax,xEnd,vStart,st,simVersion,...
        pendIndex,pos,vel,acc,theta,omega,alpha,weights);

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through experimental trials in parallel, fitting best-fit simulated
% values for each one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimizatiton function handle for calling from parfor loop
NLoptHandle = @NLopt;

% Slice start and stop indices into separate columns; this is to avoid
% making Expression1 a broadcast variable in the parfor loop
startInds = Expression1(:,1);
stopInds = Expression1(:,2);

parfor num = numStart:numEnd
    numStr = num2str(num);

    % Construct file name and load file of one experimental trial
    fileStr = strcat(subjStr,trialDate,trialStr,numStr);
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
    start = startInds(num+1);
    stop = stopInds(num+1);

    % Trim experimental data
    [pos,theta,vel,omega,acc,alpha,tdes,~] = ...
        trimData(pos,theta,vel,omega,acc,alpha,t,st,start,stop);

    % Define cart final position, initial velocity, and pendulum release
    % time from experimental data
    xEnd = pos(end);
    vStart = vel(1);
    pendIndex = find(vel > 0.1, 1);

    % Extend time vector to maximum delayed size. All time vectors will
    % be the same length, so shorter trials will be padded.
    tDesMax = tdes + delayMax;
    tc=0:st:tDesMax;

    %%% Run optimization. Fit values to selected experimental trial.
    [xopt, fmin] = feval(NLoptHandle,tc,tdes,xEnd,vStart,pendIndex,pos,vel,...
        acc,theta,omega,alpha,algorithm,setMaxEval);

    % Store best-fit parameters and objective function value in array 
    % to be accessed outside parfor-loop
    xopt_array(num+1,:) = {xopt, fmin};

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output results of optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for num = numStart:numEnd

    xopt = xopt_array{num+1,1};
    fmin = xopt_array{num+1,2};

    % Extract best-fit parameters from optimization output
    if optimizationType == "input shaping 4 impulse" || ...
            optimizationType == "input shaping 2 impulse impedance"
        optB = xopt(1);
        optK = xopt(2);
        optTcount = xopt(3);
    elseif optimizationType == "input shaping 2 impulse no impedance"
        optB = 0;
        optK = 0;
        optTcount   = xopt(1);
    end

    % Construct file name and load file of one experimental trial
    numStr = num2str(num);
    fileStr = strcat(subjStr, trialDate, trialStr, numStr);
    
    load(fileStr, 'pos', 'theta', 'vel', 'omega', 'acc', 'alpha', 't');

    % Load start and stop indices of corresponding trial. Note that trials
    % are indexed from 0, so add 1 to access correct row in array.
    start = Expression1(num+1, 1);
    stop = Expression1(num+1, 2);

    % Trim experimental data. trimData takes angles in radians and returns
    % them in degrees.
    [pos, theta, vel, omega, acc, alpha, tdes, ~] = ...
        trimData(pos, theta, vel, omega, acc, alpha, t, st, start, stop);

    % Define cart final position, initial velocity, and pendulum release
    % time from experimental data
    xEnd = pos(end);
    vStart = vel(1);
    pendIndex = find(vel > 0.1, 1);

    % Extend time vector to maximum delayed size. All time vectors will
    % be the same length, so shorter trials will be padded.
    tDesMax = tdes + delayMax;
    tc = 0:st:tDesMax;

    % Translate optimization time parameters to actual times in seconds
    optTdelay = 0.5 * (optTcount - 1) * delayMax; % Range from delayMin to delayMax
    tdessim = tdes + optTdelay;
    
    % Create system using optimal hyperparameters
    [extSys, ~, sysRigid, Td, Td1, Td2, zeta, zeta1, zeta2, ~] = ...
        sysCreate(optB, optK, intModel, impedance, printSys);

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
        [output, dur_corr, vc] = simInputShape(optB, optK, intModel, extSys,...
            sysRigid, Td1, Td2, zeta1, zeta2, tdes, tdessim, xEnd, vStart,...
            st, forwardF, simVersion, modes, pendIndex, impedance);
    elseif optimizationType == "input shaping 2 impulse impedance" || ...
            optimizationType == "input shaping 2 impulse no impedance"
        modes = 1;
        Td2 = 0;
        zeta2 = 0;
        [output, dur_corr, vc] = simInputShape(optB, optK, intModel, extSys,...
            sysRigid, Td, Td2, zeta, zeta2, tdes, tdessim, xEnd, vStart,...
            st, forwardF, simVersion, modes, pendIndex, impedance);
    end
    
    if plotPeaks
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
        % To calculate RMSE only during times when both profiles are
        % defined, shorten the longer array to the shorter array's length
        smaller = min(lensim,lenexp);
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
    dvel        = dvel / max(vel);
    squareVel   = dvel .^ 2;
    vrmse       = sqrt((1 / length(squareVel)) * sum(squareVel));
    
    % Insert trial number, fmin and optimized parameter values into array
    returnRow = num - numStart + 1;
    trial = num + 1;

    % Fill array with best-fit optimization parameters
    if optimizationType == "input shaping 4 impulse" ||...
           optimizationType == "input shaping 2 impulse impedance" 
        returnArray(returnRow,:) = [trial,fmin,vrmse,optB,optK,optTdelay];
    elseif optimizationType == "input shaping 2 impulse no impedance"
        returnArray(returnRow,:) = [trial,fmin,vrmse,optTdelay];
    end
    
    % Save best-fit output to a spreadsheet    
    durCorrVector = zeros(length(vel_sim),1);
    durCorrVector(1) = dur_corr;
    outputArray = [tc', pos_sim, theta_sim, vel_sim, omega_sim, acc_sim, alpha_sim, durCorrVector];
    profileVariables = {'Time [s]', 'Cart Pos [m]', 'Ball Angle [deg]',...
        'Cart Vel [m/s]','Ball Ang Vel [deg/s]', 'Cart Acc [m/s^2]',...
        'Ball Ang Acc [deg/s^2]', 'Pre-Trim Duration [s]'};
    ProfileT = array2table(outputArray,'VariableNames',profileVariables);
    trialNumStr = num2str(returnRow);
    fileName = strcat(folderStr, "/best-fit profiles/", trialNumStr, ".xlsx");
    writetable(ProfileT,fileName);
    
    % Check if "best" solution is actually one of filtered cases
    if mean(vel_sim)<0.1      
        disp('Chosen trajectory is 0 velocity');
    elseif max(vel_sim) > 0.5
        disp('Chosen trajectory over max velocity');
    elseif min(vel_sim) < -0.075
        disp('Chosen trajectory under min velocity');
    elseif max(theta_sim) > 50
        disp('Chosen trajectory over max ball angle');
    elseif max(omega_sim) > 75
        disp('chosen trajectory over max ball velocity');
    elseif max(acc_sim) > 5
        disp('Chosen trajectory over max cart acceleration');
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
        defaultBlue = [0, 0.4470, 0.7410];
        defaultOrange = [0.8500, 0.3250, 0.0980];

        % For characteristic trials:
        % Cart velocity (vel): plot experimental versus simulation
        set(groot,'CurrentFigure',cv);
        plot(te,vel,'LineWidth',lineWidth,'Color',defaultBlue)
        hold on;
        if plotSim
            plot(tc,vel_sim,'LineWidth',lineWidth,'Color',defaultOrange)
        end
        if plotPeaks
           plot(velChange_exp*st,velPeaks_exp,'Marker','*','MarkerSize',12,...
               'Color','magenta','LineWidth',2,'LineStyle','none')
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Time (s)','FontSize',labelFontSize)
        ylabel('Cart Velocity (m/s)','FontSize',labelFontSize)
        subjNum = erase(subjNum,'S');
        title(['Subject ', subjNum, ', Block ', block3or4 ', Trial ',...
            num2str(trialNumber)],'FontSize',titleFontSize);        
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
        end
        if sz_pos_sim(1) ~= 1
            pos_sim = pos_sim';
        end

        % Trial number as string for plot titles
        trialCaption = num2str(num+1);

        % Cart position (pos): plot experimental versus simulation
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
            if plotPeaks
                % Plot simulated local velocity max & min
                plot(velChange_sim*st,velPeaks_sim,'Marker','s',...
                    'MarkerSize',12,'Color','cyan','LineWidth',1.5,...
                    'LineStyle','none')
            end
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
        title(['Trial' ' ' trialCaption]);
        xlabel('t (s)')
        ylabel('ball acceleration (deg/s^2)')
    end    
end

% Look up NL opt algorithm name from optimization double output
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
blocknameStr = strcat("Block: ",subjStr,trialDate,trialStr);            
blocknameStr = erase(blocknameStr,"_trial_"); 
algorithm = strcat("Algorithm: ", algorithm);
Bbounds = strcat("B bounds: ",convertCharsToStrings(num2str(Bmin))," to ",...
                    convertCharsToStrings(num2str(Bmax)));
Kbounds = strcat("K bounds: ",convertCharsToStrings(num2str(Kmin))," to ",...
                    convertCharsToStrings(num2str(Kmax)));
stopTol = strcat("Stopping tolerance: ",...
    convertCharsToStrings(num2str(relative_stop)));
durationtime = strcat("Max duration error: ",...
                    convertCharsToStrings(num2str(delayMax*1000))," ms");
setMaxEval = strcat("Eval max: ",...
                    convertCharsToStrings(num2str(setMaxEval)));
wPosStr = strcat("Cart position weight: ",...
                    convertCharsToStrings(num2str(wPos)));
wThetaStr = strcat("Ball position weight: ",...
                    convertCharsToStrings(num2str(wTheta)));
wVelStr = strcat("Cart velocity weight: ",...
                    convertCharsToStrings(num2str(wVel)));
wOmegaStr = strcat("Ball velocity weight: ",...
                    convertCharsToStrings(num2str(wOmega)));
wAccStr = strcat("Cart acceleration weight: ",...
                    convertCharsToStrings(num2str(wAcc)));
wAlphaStr = strcat("Ball acceleration weight: ",...
                    convertCharsToStrings(num2str(wAlpha)));

% Add current block name and optimization settings to table
paramTable = table('Size',[numTrials,1],'VariableTypes',{'string'},...
    'VariableNames',{'Parameters'});
paramTable{1,1} = blocknameStr;
paramTable{3,1} = algorithm;
paramTable{4,1} = Bbounds;
paramTable{5,1} = Kbounds;
paramTable{6,1} = durationtime;
paramTable{7,1} = setMaxEval;
paramTable{8,1} = stopTol;
paramTable{10,1} = wPosStr;
paramTable{11,1} = wThetaStr;
paramTable{12,1} = wVelStr;
paramTable{13,1} = wOmegaStr;
paramTable{14,1} = wAccStr;
paramTable{15,1} = wAlphaStr;

if optimizationType == "input shaping 4 impulse"
    paramTable{2,1} = "4-Impulse Input Shaping Fitting";
elseif optimizationType == "input shaping 2 impulse impedance"
    paramTable{2,1} = "2-Impulse Input Shaping w/ Impedance Fitting";
elseif optimizationType == "input shaping 2 impulse no impedance"
    paramTable{2,1} = "2-Impulse Input Shaping w/out Impedance Fitting";
    % Overwrite B and K bounds with empty string
    paramTable{4,1} = " ";
    paramTable{5,1} = " ";
end 

% Convert array to table
if optimizationType == "input shaping 4 impulse" || ...
        optimizationType == "input shaping 2 impulse impedance"
    variableNames = {'Trial','Fmin','VRMSE','B','K','Tdelay'};
elseif optimizationType == "input shaping 2 impulse no impedance"
    variableNames = {'Trial','Fmin','VRMSE','Tdelay'};
end
returnTable = array2table(returnArray,'VariableNames',variableNames);

% Concatenate tables for output
T = [paramTable, returnTable];

% Analyze data in array of optimized input parameters
if optimizationType == "input shaping 4 impulse" || ...
        optimizationType == "input shaping 2 impulse impedance"
    analysisColumns = 6;
    analysisTypes = {'string','double','double','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','B','K','Tdelay'};
elseif optimizationType == "input shaping 2 impulse no impedance"
    analysisColumns = 4;
    analysisTypes = {'string','double','double','double'};
    analysisNames = {'Measure','Fmin','VRMSE','Tdelay'};
end

AnalysisT = table('Size',[size(returnTable,1), analysisColumns],...
            'VariableTypes',analysisTypes,'VariableNames',analysisNames);

% Populate analysis table with basic statistics of optimization parameters
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
fileName = strcat(folderStr, nameStr, ".xlsx");
writetable(T,fileName);
writetable(AnalysisT,fileName,'Range','I1:N6');

% Save just data without headings to a .mat file
save(strcat(topFolderStr, num2str(blockNum), ".mat"), "returnArray");

% Convert strings to character arrays for plotting
subjStr = strcat(char(subjStr),'\');
trialDate = char(trialDate);

% Give title to and save figures
if plotInd
    set(groot,'CurrentFigure',cv);
    saveas(cv, strcat(folderStr, "cart_vel_", nameStr, ".png"));
else
    figTitle = ['Experimental vs Best Fit Cart Position Profiles,' ' '...
        subjStr trialDate];
    sgtitle(cp,figTitle,'FontSize',16)
    
    figTitle = ['Experimental vs Best Fit Ball Angle Profiles,' ' '...
        subjStr trialDate];
    sgtitle(bp,figTitle,'FontSize',16)

    figTitle = ['Experimental vs Best Fit Cart Velocity Profiles: Subject ',...
        subjNum, ' Block ', block3or4 ];
    sgtitle(cv,figTitle,'FontSize',16)

    figTitle = ['Experimental vs Best Fit Ball Angular Velocity Profiles,'...
        ' ' subjStr trialDate];
    sgtitle(bv,figTitle,'FontSize',16)

    figTitle = ['Experimental vs Best Fit Cart Acceleration Profiles,' ' '...
        subjStr trialDate];
    sgtitle(ca,figTitle,'FontSize',16)

    figTitle = ['Experimental vs Best Fit Ball Angular Acceleration Profiles,'...
        ' ' subjStr trialDate];
    sgtitle(ba,figTitle,'FontSize',16)

    saveas(cp, strcat(folderStr, "cart_pos_", nameStr, ".png"))      
    saveas(bp, strcat(folderStr, "ball_angle_", nameStr, ".png"))      
    saveas(cv, strcat(folderStr, "cart_vel_", nameStr, ".png"))
    saveas(bv, strcat(folderStr, "ball_angular_vel_", nameStr, ".png"))
    saveas(ca, strcat(folderStr, "cart_acc_", nameStr, ".png"))
    saveas(ba, strcat(folderStr, "ball_angular_acc_", nameStr, ".png"))
end

% Give focus back to command window when saving is complete
commandwindow

end