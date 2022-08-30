% minVelStats

load('best fit simulation/correlationData.mat')

% Extract correlation coefficients for different simulations and
% experimental data into a matrix
rMatrix = cell2mat(rBySubject(2:12,2:16));
rExp = rMatrix(:,15);
simLabels = rBySubject(1,2:15);

pairedTestResults = { 'Data Type';...
                'Significant to 5%?';...
                '5% p value';...
                'Significant to 2%?';...
                '2% p value'};
            
singleTestResults = pairedTestResults;

disp('Paired t-test results')

% Conduct paired t-test to determine if mean of correlation coefficients
% are significantly different from each other at 5% and 2% signifiance levels
for simNum = 1:14
    simLabel = rBySubject{1,simNum+1}
    
    tail = 'both';
    [h5,p5,ci,stats] = ttest(rExp,rMatrix(:,simNum),'Tail',tail)

    [h2,p2,ci,stats] = ttest(rExp,rMatrix(:,simNum),'Tail',tail,'Alpha',.02)
    
    pairedTestResults(1,simNum+1) = {simLabels{simNum}};
    pairedTestResults(2:5,simNum+1) = {h5; p5; h2; p2};
    
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Single t-test results')

% Conduct single t-test to determine if mean of correlation coefficients
% are significantly less than 0 at 5% and 2% signifiance levels
for simNum = 1:14
    simLabel = rBySubject{1,simNum+1}
    
    tail = 'left';
    [h5,p5,ci,stats] = ttest(rMatrix(:,simNum),0,'Tail',tail)

    [h2,p2,ci,stats] = ttest(rMatrix(:,simNum),0,'Tail',tail)
    
    singleTestResults(1,simNum+1) = {simLabels{simNum}};
    singleTestResults(2:5,simNum+1) = {h5; p5; h2; p2};
    
end
