function [subjNum, subjStr, trialDate, trialStr, blockStr, plotColor,...
        minVelOutliers, invalidTrials] = blockDictionary(blockNum)
% BLOCKDICTIONARY
%
% Contains information about all experimental blocks as well as unique
% colors to plot for each subject. Can be used by multiple functions to
% loop through all experimental blocks.
%
% The minVelOutliers array contains the trial numbers that are outliers
% when looking at the correlation of inter-peak minimum velocity and trial
% duration.

            if blockNum == 1
                % S1: EF Block 3
                subjNum     = 'S1';
                subjStr     = 'EF';
                trialDate   = '_09Jul2019';
                trialStr    = '_14-20-29_trial_';
                blockStr    = 'EFb3.mat';
                plotColor   = [0.368417 0.506779 0.709798];
                minVelOutliers = [1, 2, 6, 8, 9, 10, 15, 16, 17, 18, 19,...
                    20, 23, 26, 27, 31, 45, 46];
                invalidTrials = [];

            elseif blockNum == 2
                % S1: EF Block 4
                subjNum     = 'S1';
                subjStr     = 'EF';
                trialDate   = '_09Jul2019';
                trialStr    = '_14-30-08_trial_';
                blockStr    = 'EFb4.mat';
                plotColor   = [0.368417 0.506779 0.709798];
                minVelOutliers = [9, 12, 13, 17, 19, 21, 22, 27, 29, 33,...
                    35, 40];
                invalidTrials = [];
    
            elseif blockNum == 3
                % S2: NO Block 3
                subjNum     = 'S2';
                subjStr     = 'NO';
                trialDate   = '_15Jul2019';
                trialStr    = '_13-17-36_trial_';
                blockStr    = 'NOb3.mat';
                plotColor   = [0.880722 0.611041 0.142051];
                minVelOutliers = [];
                invalidTrials = [];

            elseif blockNum == 4
                % S2: NO Block 4
                subjNum     = 'S2';
                subjStr     = 'NO';
                trialDate   = '_15Jul2019';
                trialStr    = '_13-27-13_trial_';
                blockStr    = 'NOb4.mat';
                plotColor   = [0.880722 0.611041 0.142051];
                minVelOutliers = [4, 50];
                invalidTrials = [];

            elseif blockNum == 5
                % S3: MHH Block 3
                subjNum     = 'S3';
                subjStr     = 'MH';
                trialDate   = '_16Jul2019';
                trialStr    = '_10-04-29_trial_';
                blockStr    = 'MHHb3.mat';
                plotColor   = [0.560181 0.691569 0.194885];
                minVelOutliers = [3, 15, 42, 46, 47, 50];
                invalidTrials = [];

            elseif blockNum == 6
                % S3: MHH Block 4
                subjNum     = 'S3';
                subjStr     = 'MH';
                trialDate   = '_16Jul2019';
                trialStr    = '_10-14-05_trial_';
                blockStr    = 'MHHb4.mat';
                plotColor   = [0.560181 0.691569 0.194885];
                minVelOutliers = [1, 3, 6, 7, 8, 11, 16, 33, 41];
                invalidTrials = [];

            elseif blockNum == 7
                % S4: CM Block 3
                subjNum     = 'S4';
                subjStr     = 'CM';
                trialDate   = '_17Jul2019';
                trialStr    = '_16-37-24_trial_';
                blockStr    = 'CMb3.mat';
                plotColor   = [0.922526,0.385626,0.209179];
                minVelOutliers = [1, 6, 31, 33, 38];
                invalidTrials = [33];
                

            elseif blockNum == 8
                % S4: CM Block 4
                subjNum     = 'S4';
                subjStr     = 'CM';
                trialDate   = '_17Jul2019';
                trialStr    = '_16-47-55_trial_';
                blockStr    = 'CMb4.mat';
                plotColor   = [0.922526,0.385626,0.209179];
                minVelOutliers = [42];
                invalidTrials = [];

            elseif blockNum == 9
                % S5: AP Block 3
                subjNum     = 'S5';
                subjStr     = 'AP';
                trialDate   = '_19Jul2019';
                trialStr    = '_13-31-43_trial_';
                blockStr    = 'APb3.mat';
                plotColor   = [0.528488,0.470624,0.701351];
                minVelOutliers = [16];
                invalidTrials = [];

            elseif blockNum == 10
                % S5: AP Block 4
                subjNum     = 'S5';
                subjStr     = 'AP';
                trialDate   = '_19Jul2019';
                trialStr    = '_13-41-46_trial_';
                blockStr    = 'APb4.mat';
                plotColor   = [0.528488,0.470624,0.701351];
                minVelOutliers = [];
                invalidTrials = [];

            elseif blockNum == 11
                % S6: DD Block 3
                subjNum     = 'S6';
                subjStr     = 'DD';
                trialDate   = '_22Jul2019';
                trialStr    = '_09-15-05_trial_';
                blockStr    = 'DDb3.mat';
                plotColor   = [0.772079,0.431554,0.102387];
                minVelOutliers = [6, 7, 8, 9, 11, 12, 25];
                invalidTrials = [];

            elseif blockNum == 12
                % S6: DD Block 4
                subjNum     = 'S6';
                subjStr     = 'DD';
                trialDate   = '_22Jul2019';
                trialStr    = '_09-25-00_trial_';
                blockStr    = 'DDb4.mat';
                plotColor   = [0.772079,0.431554,0.102387];
                minVelOutliers = [6, 7, 47];
                invalidTrials = [];

            elseif blockNum == 13
                % S7: SC Block 3
                subjNum     = 'S7';
                subjStr     = 'SC';
                trialDate   = '_23Jul2019';
                trialStr    = '_11-19-32_trial_';
                blockStr    = 'SCb3.mat';
                plotColor   = [0.363898,0.618501,0.782349];
                minVelOutliers = [2, 14, 32];
                invalidTrials = [];

            elseif blockNum == 14
                % S7: SC Block 4
                subjNum     = 'S7';
                subjStr     = 'SC';
                trialDate   = '_23Jul2019';
                trialStr    = '_11-29-07_trial_';
                blockStr    = 'SCb4.mat';
                plotColor   = [0.363898,0.618501,0.782349];
                minVelOutliers = [23, 33];
                invalidTrials = [];
            
            elseif blockNum == 15
                % S8: AK Block 3
                subjNum     = 'S8';
                subjStr     = 'AK';
                trialDate   = '_11Nov2020';
                trialStr    = '_15-27-54_trial_';
                blockStr    = 'AKb3.mat';
                plotColor   = [1,0.75,0];
                minVelOutliers = [1, 4, 16, 19, 29, 50];
                invalidTrials = [];
            
            elseif blockNum == 16
                % S8: AK Block 4
                subjNum     = 'S8';
                subjStr     = 'AK';
                trialDate   = '_11Nov2020';
                trialStr    = '_15-38-17_trial_';
                blockStr    = 'AKb4.mat';
                plotColor   = [1,0.75,0];
                minVelOutliers = [1, 19, 22, 45];
                invalidTrials = [];
            
            elseif blockNum == 17
                % S9: GB Block 3
                subjNum     = 'S9';
                subjStr     = 'GB';
                trialDate   = '_12Nov2020';
                trialStr    = '_14-19-04_trial_';
                blockStr    = 'GBb3.mat';
                plotColor   = [0.647624,0.37816,0.614037];
                minVelOutliers = [1, 27, 43, 48];
                invalidTrials = [];
            
            elseif blockNum == 18
                % S9: GB Block 4
                subjNum     = 'S9';
                subjStr     = 'GB';
                trialDate   = '_12Nov2020';
                trialStr    = '_14-29-34_trial_';
                blockStr    = 'GBb4.mat';
                plotColor   = [0.647624,0.37816,0.614037];
                minVelOutliers = [8, 39, 50];
                invalidTrials = [];
            
            elseif blockNum == 19
                % S10: HC Block 3
                subjNum     = 'S10';
                subjStr     = 'HC';
                trialDate   = '_13Nov2020';
                trialStr    = '_16-32-53_trial_';
                blockStr    = 'HCb3.mat';
                plotColor   = [0.571589,0.586483,0];
                minVelOutliers = [24, 48];
                invalidTrials = [];
             
            elseif blockNum == 20
                % S10: HC Block 4
                subjNum     = 'S10';
                subjStr     = 'HC';
                trialDate   = '_13Nov2020';
                trialStr    = '_16-43-21_trial_';
                blockStr    = 'HCb4.mat';
                plotColor   = [0.571589,0.586483,0];
                minVelOutliers = [2, 45, 50];
                invalidTrials = [2];
            
            elseif blockNum == 21
                % S11: JM Block 3
                subjNum     = 'S11';
                subjStr     = 'JM';
                trialDate   = '_13Nov2020';
                trialStr    = '_13-21-05_trial_';
                blockStr    = 'JMb3.mat';
                plotColor   = [0.915,0.3325,0.2125];
                minVelOutliers = [1, 5, 12, 31, 37, 48];
                invalidTrials = [];
                
            elseif blockNum == 22
                % S11: JM Block 4
                subjNum     = 'S11';
                subjStr     = 'JM';
                trialDate   = '_13Nov2020';
                trialStr    = '_13-31-11_trial_';
                blockStr    = 'JMb4.mat';
                plotColor   = [0.915,0.3325,0.2125];
                minVelOutliers = [38, 45, 46];
                invalidTrials = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Below this line is information for Block 2 %
            elseif blockNum == 23
                % S1: EF Block 2
                subjNum     = 'S1';
                subjStr     = 'EF';
                trialDate   = '_09Jul2019';
                trialStr    = '_14-10-49_trial_';
                blockStr    = 'EFb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
                
            elseif blockNum == 24
                % S2: NO Block 2
                subjNum     = 'S2';
                subjStr     = 'NO';
                trialDate   = '_15Jul2019';
                trialStr    = '_13-07-18_trial_';
                blockStr    = 'NOb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
                
            elseif blockNum == 25
                % S3: MHH Block 2
                subjNum     = 'S3';
                subjStr     = 'MH';
                trialDate   = '_15Jul2019';
                trialStr    = '_14-05-11_trial_';
                blockStr    = 'MHHb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [20];
                
            elseif blockNum == 26
                % S4: CM Block 2
                subjNum     = 'S4';
                subjStr     = 'CM';
                trialDate   = '_17Jul2019';
                trialStr    = '_16-27-10_trial_';
                blockStr    = 'CMb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [1];
                
            elseif blockNum == 27
                % S5: AP Block 2
                subjNum     = 'S5';
                subjStr     = 'AP';
                trialDate   = '_19Jul2019';
                trialStr    = '_13-22-10_trial_';
                blockStr    = 'APb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
                
            elseif blockNum == 28
                % S6: DD Block 2
                subjNum     = 'S6';
                subjStr     = 'DD';
                trialDate   = '_22Jul2019';
                trialStr    = '_09-06-17_trial_';
                blockStr    = 'DDb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
                
            elseif blockNum == 29
                % S7: SC Block 2
                subjNum     = 'S7';
                subjStr     = 'SC';
                trialDate   = '_23Jul2019';
                trialStr    = '_11-05-44_trial_';
                blockStr    = 'SCb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [35];
                
            elseif blockNum == 30
                % S8: AK Block 2
                subjNum     = 'S8';
                subjStr     = 'AK';
                trialDate   = '_11Nov2020';
                trialStr    = '_15-17-57_trial_';
                blockStr    = 'AKb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [5,10];
                
            elseif blockNum == 31
                % S9: GB Block 2
                subjNum     = 'S9';
                subjStr     = 'GB';
                trialDate   = '_12Nov2020';
                trialStr    = '_14-08-56_trial_';
                blockStr    = 'GBb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [1, 2];
                
            elseif blockNum == 32
                % S10: HC Block 2
                subjNum     = 'S10';
                subjStr     = 'HC';
                trialDate   = '_13Nov2020';
                trialStr    = '_16-21-56_trial_';
                blockStr    = 'HCb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
                
            elseif blockNum == 33
                % S11: JM Block 2
                subjNum     = 'S11';
                subjStr     = 'JM';
                trialDate   = '_13Nov2020';
                trialStr    = '_13-12-10_trial_';
                blockStr    = 'JMb2.mat';
                plotColor   = 	[0 0 0];
                minVelOutliers = [];
                invalidTrials = [];
            end            
end