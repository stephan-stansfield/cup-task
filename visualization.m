% figTrialOverlayVisualization
%
% Plots experimental and simulated trials from one subject in a single 
% overlay plot. Can plot individual trials from one or both test blocks,
% and the average of all plotted trials.

function visualization()
    
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

% Choose whether to plot experimental trials (true) or just
% simulated trials (false)
plotExp = true;

% Choose whether to plot just experimental trials (true)
plotExpOnly = false;

% Choose whether to plot individual simulated trials (true)
plotIndividualSim = false;

% Choose whether to plot both blocks 3 & 4 for each subject (true) or
% to plot each block individually (false)
pairs = true;

% Choose whether to save generated figures
save = false;

% Choose whether to skip trials without 2 clear velocity peaks (true)
% or to visualize all valid trials (false)
peaksOnly = false;

% Choose whether to plot figures with subject number as title (true)
% or with no title (false)
subjTitle = true;

% Define colors for plots
startBlack  = [0.5, 0.5, 0.5];
startRed    = [255, 153, 153]/255;
startOrange = [255, 194, 153]/255;
startPink   = [255, 153, 255]/255;
startGreen  = [204, 255, 204]/255;
startBlue   = [153, 194, 255]/255;
startPurple = [204, 153, 255]/255;

% Was using bolder colors for experimental average but changed to black.
% Keep values here for handy reference if want to change back.
% boldRed     = [255, 0, 0]/255;
% boldOrange  = [255, 102, 0]/255;
% boldPink    = [255, 0, 255]/255;
% boldGreen   = [31, 122, 31]/255;
% boldBlue    = [0, 102, 255]/255;
% boldPurple  = [191, 0, 255]/255;

% Choose which experimental blocks to visualize (all acceptable blocks
% would be 1-22). See "blockDictionary.m" for a key to which blocks
% correspond with which subject number.
blockStart  = 17;
blockEnd    = 18;
numBlocks   = blockEnd - blockStart + 1;

% Choose range of experimental trials to visualize
numStart    = 1;
numEnd      = 50;

% Cycle through simulation types
for count = 13
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

    for blockNum = blockStart:blockEnd

        % Get information about experimental block
        [subjNum, subjStr, trialDate, trialStr, blockStr, ~,~,...
            invalidTrials] = blockDictionary(blockNum);
        
        block3 = mod(blockNum,2);
        block4 = ~mod(blockNum,2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Start visualization function
        
        load(blockStr,'Expression1')                                                % Load trimming indices
        
        if pairs
            numTrials = 100;
        else
            numTrials = 50;
        end
        
        if pairs && block4
            % If plotting both blocks and on block 4, don't initialize
        else
            % If plotting both blocks and on block 3, or if plotting
            % individual blocks, initialize arrays to hold all data
            posArray    = zeros(numTrials,normLen);
            thetaArray  = zeros(numTrials,normLen);
            velArray    = zeros(numTrials,normLen);
            omegaArray  = zeros(numTrials,normLen);
            accArray    = zeros(numTrials,normLen);
            alphaArray  = zeros(numTrials,normLen);
            tdesArray   = zeros(numTrials,1);

            posSimArray    = zeros(numTrials,normLen);
            thetaSimArray  = zeros(numTrials,normLen);
            velSimArray    = zeros(numTrials,normLen);
            omegaSimArray  = zeros(numTrials,normLen);
            accSimArray    = zeros(numTrials,normLen);
            alphaSimArray  = zeros(numTrials,normLen);
            tdesSimArray   = zeros(numTrials,1);
        end

        % Folder holding best-fit trials from current block
        simFolder = strcat("best fit simulation/",parentFolder,subjStr,...
            trialDate,trialStr);
        simFolder = erase(simFolder, "_trial_");

        if pairs && block4
            % don't initialize figure windows
        else
            % Initialize figure windows
            figure('Name','Cart position','units','normalized');
            cp = get(gcf,'Number');

            figure('Name','Cart velocity','units','normalized');
            cv = get(gcf,'Number');

            figure('Name','Cart acceleration','units','normalized');
            ca = get(gcf,'Number');

            figure('Name','Ball angle','units','normalized');
            bp = get(gcf,'Number');

            figure('Name','Ball angular velocity','units','normalized');
            bv = get(gcf,'Number');

            figure('Name','Ball angular acceleration','units','normalized');
            ba = get(gcf,'Number');
        end

        % Loop through trials, plotting each one and storing values in array
        for num = numStart:1:numEnd
            
            if peaksOnly
                % Skip invalid trials & trials that didn't have 2 clear peaks
                if ismember(num,minVelOutliers)
                    continue
                end
            end
            
            % Skip invalid trials
            if ismember(num,invalidTrials)
                continue
            end
            
            numStr = num2str(num-1);

            % Construct file name and load file of one experimental trial
            fileStr = strcat(subjStr,trialDate,trialStr,numStr);
            load(fileStr,'pos','theta','vel','omega','acc','alpha','t')

            % Load start and stop indices of corresponding trial. Note that trials
            % are indexed from 0, so add 1 to access correct row in array.
            start   = Expression1(num,1);
            stop    = Expression1(num,2);

            % Trim experimental data
            [pos,theta,vel,omega,acc,alpha,tdes] = trimData(pos,theta,vel,...
                omega,acc,alpha,t,st,start,stop);

            % Normalize experimental profiles
            [pos,theta,vel,omega,acc,alpha] = normalize(pos,theta,vel,...
                omega,acc,alpha,tdes,st,stNorm,[]);

            % Load best-fit SIMULATED trial
            simProfiles = readmatrix(strcat(simFolder, "/best-fit profiles/", num2str(num), ".xlsx"));
            tdesSim     = simProfiles(end,1);
            posSim      = simProfiles(:,2);
            thetaSim    = simProfiles(:,3);
            velSim      = simProfiles(:,4);
            omegaSim    = simProfiles(:,5);
            accSim      = simProfiles(:,6);
            alphaSim    = simProfiles(:,7);

            % Normalize simulated profiles
            [posSim,thetaSim,velSim,omegaSim,accSim,alphaSim] = ...
                normalize(posSim,thetaSim,velSim,omegaSim,accSim,alphaSim,...
                tdesSim,st,stNorm,[]);

            if pairs
                if block4
                    arrayNum = num+50;
                else
                    arrayNum = num;
                end
            else
                arrayNum = num;
            end
            
            % Add normalized profiles to arrays with all trials from block
            posArray(arrayNum,:)      = pos;
            thetaArray(arrayNum,:)    = theta;
            velArray(arrayNum,:)      = vel;
            omegaArray(arrayNum,:)    = omega;
            accArray(arrayNum,:)      = acc;
            alphaArray(arrayNum,:)    = alpha;
            tdesArray(arrayNum,:)     = tdes;

            posSimArray(arrayNum,:)   = posSim;
            thetaSimArray(arrayNum,:) = thetaSim;
            velSimArray(arrayNum,:)   = velSim;
            omegaSimArray(arrayNum,:) = omegaSim;
            accSimArray(arrayNum,:)   = accSim;
            alphaSimArray(arrayNum,:) = alphaSim;
            tdesSimArray(arrayNum,:)  = tdesSim;

            if plotExp
                %%%% Add normalized experimental trials to plot %%%%
                % Multiplier for color gradient
                if pairs
                    multiplier = 0.75+0.25*arrayNum/100;  
                else
                    multiplier = 0.75+0.25*arrayNum/50;
                end

                % Cart position
                set(groot,'CurrentFigure',cp);
                plot(tNorm,pos,'Color',startGreen*multiplier);
                hold on;

                % Cart velocity
                set(groot,'CurrentFigure',cv);
                plot(tNorm,vel,'Color',startBlue*multiplier);
                hold on;

                % Cart acceleration
                set(groot,'CurrentFigure',ca);
                plot(tNorm,acc,'Color',startPurple*multiplier);
                hold on;

                % Ball angle
                set(groot,'CurrentFigure',bp);
                plot(tNorm,theta,'Color',startRed*multiplier);
                hold on;

                % Ball velocity
                set(groot,'CurrentFigure',bv);
                plot(tNorm,omega,'Color',startOrange*multiplier);
                hold on;

                % Ball acceleration
                set(groot,'CurrentFigure',ba);
                plot(tNorm,alpha,'Color',startPink*multiplier);
                hold on;
            end
            
            if plotIndividualSim
                %%%% Add normalized simulated trials to plot %%%%
                % Multiplier for color gradient
                if pairs
                    multiplier = 0.75+0.25*arrayNum/100;  
                else
                    multiplier = 0.75+0.25*arrayNum/50;
                end

                % Cart velocity
                set(groot,'CurrentFigure',cv);
                plot(tNorm,velSim,'Color',startBlue*multiplier);
                hold on;
            end

        end
        
        % If plotting both blocks and currently on block 3, continue onto
        % block 4 before taking average of all trials
        if pairs
            if block3
                continue
            end
        end
        
        % Delete skipped trials from arrays. Filter them out by finding the
        % trials whose final position is 0.
        [zeroRow, ~] = find(~posSimArray(:,end));
        posArray(zeroRow,:)     = [];
        thetaArray(zeroRow,:)   = [];
        velArray(zeroRow,:)     = [];
        omegaArray(zeroRow,:)   = [];
        accArray(zeroRow,:)     = [];
        alphaArray(zeroRow,:)   = [];
        
        posSimArray(zeroRow,:)     = [];
        thetaSimArray(zeroRow,:)   = [];
        velSimArray(zeroRow,:)     = [];
        omegaSimArray(zeroRow,:)   = [];
        accSimArray(zeroRow,:)     = [];
        alphaSimArray(zeroRow,:)   = [];

        %%%% Plot average of experimental and simulated trials %%%%
        % Take average of all trials
        posAvg      = mean(posArray);
        thetaAvg    = mean(thetaArray);
        velAvg      = mean(velArray);
        omegaAvg    = mean(omegaArray);
        accAvg      = mean(accArray);
        alphaAvg    = mean(alphaArray);

        posSimAvg   = mean(posSimArray);
        thetaSimAvg = mean(thetaSimArray);
        velSimAvg   = mean(velSimArray);
        omegaSimAvg = mean(omegaSimArray);
        accSimAvg   = mean(accSimArray);
        alphaSimAvg = mean(alphaSimArray);

        % Set plot properties
        % Experimental trials average
        expWidth = 3;
        expStyle = '--';

        % Simulated trials average
        simWidth = 3;
        simStyle = '-';

        % Figure properties
        labelFontSize = 20;
        axisFontSize = 20;
        titleFontSize = 28;
        if subjTitle
            titleStr = subjNum;
            if plotIndividualSim
%                 titleStr = [subjNum, ': Individual Simulated Trials and Average'];
                titleFontSize = 20;
            end
        else
            titleStr = "";
        end

        % Cart position
        set(groot,'CurrentFigure',cp);
        if plotExp
            plot(tNorm,posAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,posSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Cart Position (m)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        % Cart velocity
        set(groot,'CurrentFigure',cv);
        if plotExp
            plot(tNorm,velAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,velSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Cart Velocity (m/s)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        % Cart acceleration
        set(groot,'CurrentFigure',ca);
        if plotExp
            plot(tNorm,accAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,accSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Cart Acceleration (m/s^2)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        % Ball angle
        set(groot,'CurrentFigure',bp);
        if plotExp
            plot(tNorm,thetaAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,thetaSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Ball Angle (deg)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        % Ball velocity
        set(groot,'CurrentFigure',bv);
        if plotExp
            plot(tNorm,omegaAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,omegaSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Ball Angular Velocity (deg/s)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        % Ball acceleration
        set(groot,'CurrentFigure',ba);
        if plotExp
            plot(tNorm,alphaAvg,'Color','black','LineWidth',expWidth,'LineStyle',expStyle);
            hold on;
        end
        if ~plotExpOnly
            plot(tNorm,alphaSimAvg,'Color','black','LineWidth',simWidth,'LineStyle',simStyle);
        end
        set(gca,'FontSize',axisFontSize)
        xlabel('Normalized time','FontSize',labelFontSize)
        ylabel('Ball Angular Acceleration (deg/s^2)','FontSize',labelFontSize)
        title(titleStr,'FontSize',titleFontSize)
        hold off;

        if save
            % Save figures as .png files
            cpStr       = "1 cart_pos_";
            bpStr       = "2 ball_angle_";
            cvStr       = "3 cart_vel_";
            bvStr       = "4 ball_angular_vel_";
            caStr       = "5 cart_acc_";
            baStr       = "6 ball_angular_acc_";

            if subjTitle
                nameStr = strcat(subjNum,". ",subjStr,trialDate,trialStr);
            else
                nameStr = strcat(subjNum,". ",subjStr,trialDate,trialStr,"_No Titles");
            end
            nameStr = erase(nameStr, "_trial_");
            
            % Make visualization output folder within folder holding
            % best-fit simulation trials of one type
            if(not(isfolder(strcat("best fit simulation/",parentFolder,"_Visualization Output"))))
                mkdir(strcat("best fit simulation/",parentFolder),"_Visualization Output");
            end
            
            % Make folder to hold experimental trial-only visualizations
            if(not(isfolder(strcat("best fit simulation/","_Experimental Visualization"))))
                mkdir(strcat("best fit simulation/"),"_Experimental Visualization");
            end
            
            if peaksOnly
               nameStr = strcat(nameStr, "_Peaks Only"); 
            end
            
            if plotExpOnly
                mkdir(strcat("best fit simulation/","_Experimental Visualization"),nameStr);
                saveFolder = strcat("best fit simulation/","_Experimental Visualization/",nameStr,"/");
            elseif plotIndividualSim && ~plotExp
                mkdir(strcat("best fit simulation/","_Simulated Only Visualization"),nameStr);
                saveFolder = strcat("best fit simulation/","_Simulated Only Visualization/",nameStr,"/");
            else
                mkdir(strcat("best fit simulation/",parentFolder,"_Visualization Output"),nameStr);
                saveFolder = strcat("best fit simulation/",parentFolder,"_Visualization Output/",nameStr,"/");
            end

            saveas(cp, strcat(saveFolder, cpStr, nameStr, ".png"))          
            saveas(bp, strcat(saveFolder, bpStr, nameStr, ".png"))                                                            
            saveas(cv, strcat(saveFolder, cvStr, nameStr, ".png"))
            saveas(bv, strcat(saveFolder, bvStr, nameStr, ".png"))
            saveas(ca, strcat(saveFolder, caStr, nameStr, ".png"))
            saveas(ba, strcat(saveFolder, baStr, nameStr, ".png"))

            close all
        end
        
        
        % End visualization function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
    
end
    
end
        