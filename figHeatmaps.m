% figHeatmaps

clear all;
close all;

saveFile = true;

% Define labels
%% NOTE! CHANGE ORDER BACK IF EDITING HEATMAPS (OR GO THROUGH AND DO DATA IN THIS NEW ORDER BLEOW %%
Simulation = {'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'; ...
    'IS'; 'MM-NF'; 'NI-NF'; 'RB-NF'; 'SM-NF'; 'FM-NF'; 'MM-FF';...
    'NI-FF'; 'RB-FF'; 'SM-FF'; 'FM-FF';'SB'; 'SB-NF'; 'SB-FF'};

Subject = {'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1'; 'S1';...
    'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; 'S2'; ...
    'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; 'S3'; ...
    'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; 'S4'; ...
    'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; 'S5'; ...
    'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; 'S6'; ...
    'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; 'S7'; ...
    'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; 'S8'; ...
    'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; 'S9'; ...
    'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; 'S10'; ...
    'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'; 'S11'};

% Load saved VAF and impedance data
load('best fit simulation/VAFarrays.mat')
load('best fit simulation/best fit impedance.mat')

% Unroll matrices into vectors
meanVAFvector = meanVAF(:);
stdVAFvector = sdevVAF(:);
meanKvector = meanK(:);
stdKvector = stdK(:);
meanBvector = meanB(:);
stdBvector = stdB(:);

% Convert VAF double values to percentages
% put code here later

% Convert to tables
meanVAFtable = table(Simulation,Subject,meanVAFvector);
stdVAFtable = table(Simulation,Subject,stdVAFvector);
meanKtable = table(Simulation,Subject,meanKvector);
stdKtable = table(Simulation,Subject,stdKvector);
meanBtable = table(Simulation,Subject,meanBvector);
stdBtable = table(Simulation,Subject,stdBvector);

% Create heatmaps
meanVAFfig = figure();
hMeanVAF = heatmap(meanVAFtable,'Subject','Simulation','ColorVariable','meanVAFvector',...
    'ColorMethod','none','Colormap',hot);
stdVAFfig = figure();
hStdVAF = heatmap(stdVAFtable,'Subject','Simulation','ColorVariable','stdVAFvector',...
    'ColorMethod','none','Colormap',hot);

meanKfig = figure();
hMeanK = heatmap(meanKtable,'Subject','Simulation','ColorVariable','meanKvector',...
    'ColorMethod','none','Colormap',copper);
stdKfig = figure();
hStdK = heatmap(stdKtable,'Subject','Simulation','ColorVariable','stdKvector',...
    'ColorMethod','none','Colormap',copper);

meanBfig = figure();
hMeanB = heatmap(meanBtable,'Subject','Simulation','ColorVariable','meanBvector',...
    'ColorMethod','none','Colormap',summer);
stdBfig = figure();
hStdB = heatmap(stdBtable,'Subject','Simulation','ColorVariable','stdBvector',...
    'ColorMethod','none','Colormap',summer);

% Re-order subjects and simulation types (default is alphabetical) and
% re-title heatmaps
subjectOrder = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};
simulationOrder = {'IS', 'MM-FF', 'NI-FF', 'RB-FF', 'SM-FF', 'FM-FF',...
    'MM-NF', 'NI-NF', 'RB-NF', 'SM-NF', 'FM-NF', 'SB', 'SB-FF', 'SB-NF'};

% VAF
hMeanVAF.SourceTable.Subject = categorical(hMeanVAF.SourceTable.Subject);
hMeanVAF.SourceTable.Subject = reordercats(hMeanVAF.SourceTable.Subject,subjectOrder);
hMeanVAF.SourceTable.Simulation = categorical(hMeanVAF.SourceTable.Simulation);
hMeanVAF.SourceTable.Simulation = reordercats(hMeanVAF.SourceTable.Simulation,simulationOrder);
hMeanVAF.Title = 'Mean Velocity VAF by Simulation Type and Subject';

hStdVAF.SourceTable.Subject = categorical(hStdVAF.SourceTable.Subject);
hStdVAF.SourceTable.Subject = reordercats(hStdVAF.SourceTable.Subject,subjectOrder);
hStdVAF.SourceTable.Simulation = categorical(hStdVAF.SourceTable.Simulation);
hStdVAF.SourceTable.Simulation = reordercats(hStdVAF.SourceTable.Simulation,simulationOrder);
hStdVAF.Title = 'Velocity VAF Standard Deviation by Simulation Type and Subject';

% K
hMeanK.SourceTable.Subject = categorical(hMeanK.SourceTable.Subject);
hMeanK.SourceTable.Subject = reordercats(hMeanK.SourceTable.Subject,subjectOrder);
hMeanK.SourceTable.Simulation = categorical(hMeanK.SourceTable.Simulation);
hMeanK.SourceTable.Simulation = reordercats(hMeanK.SourceTable.Simulation,simulationOrder);
hMeanK.Title = 'Mean Hand Stiffness by Simulation Type and Subject';

hStdK.SourceTable.Subject = categorical(hStdK.SourceTable.Subject);
hStdK.SourceTable.Subject = reordercats(hStdK.SourceTable.Subject,subjectOrder);
hStdK.SourceTable.Simulation = categorical(hStdK.SourceTable.Simulation);
hStdK.SourceTable.Simulation = reordercats(hStdK.SourceTable.Simulation,simulationOrder);
hStdK.Title = 'Hand Stiffness Standard Deviation by Simulation Type and Subject';

% B
hMeanB.SourceTable.Subject = categorical(hMeanB.SourceTable.Subject);
hMeanB.SourceTable.Subject = reordercats(hMeanB.SourceTable.Subject,subjectOrder);
hMeanB.SourceTable.Simulation = categorical(hMeanB.SourceTable.Simulation);
hMeanB.SourceTable.Simulation = reordercats(hMeanB.SourceTable.Simulation,simulationOrder);
hMeanB.Title = 'Mean Hand Damping by Simulation Type and Subject';

hStdB.SourceTable.Subject = categorical(hStdB.SourceTable.Subject);
hStdB.SourceTable.Subject = reordercats(hStdB.SourceTable.Subject,subjectOrder);
hStdB.SourceTable.Simulation = categorical(hStdB.SourceTable.Simulation);
hStdB.SourceTable.Simulation = reordercats(hStdB.SourceTable.Simulation,simulationOrder);
hStdB.Title = 'Hand Damping Standard Deviation by Simulation Type and Subject';

%%%%
if saveFile
   Folder = 'best fit simulation/_figures/';
    
   saveas(meanVAFfig, [Folder, 'mean VAF.png']);
   saveas(stdVAFfig, [Folder, 'std dev VAF.png']);
   saveas(meanKfig, [Folder, 'mean K.png']);
   saveas(stdKfig, [Folder, 'std dev K.png']);
   saveas(meanBfig, [Folder, 'mean B.png']);
   saveas(stdBfig, [Folder, 'std dev B.png']);
    
end


