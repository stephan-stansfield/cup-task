% Calculate statistical significance of difference between model mean
% values using a t-test.

clear all
close all
clc

% %% Looking at mean of subject **median** VAFs
% C = readcell('data analysis/VAF summary 2023-03-26.xlsx'); % simulations for 2023 journal paper
% 
% % Read all columns:
% x = C(1,2:end);                     % x-axis labels
% u = 100*cell2mat(C(2,2:end));       % mean VAFs
% s = 100*cell2mat(C(3,2:end));       % standard deviations

% load array containing mean VAF by simulation (rows) and subject (columns)
load('statistics/statisticsVAF.mat');

modelNum = 5;               % number of models to compare against original 
pArray = ones(1, modelNum); % array to hold all p values
h5Array = ones(1, modelNum); % array to hold true/false results to 5% level
h1Array = ones(1, modelNum); % array to hold true/false results to 1% level

% compare each control model to original input shaping
for model = 1:5
    [h5, p] = ttest(meanVAF(1, :), meanVAF(model+1, :));
    pArray(model) = p;
    h5Array(model) = h5;

    [h1, p] = ttest(meanVAF(1, :), meanVAF(model+1, :), 'Alpha', 0.01);
    h1Array(model) = h1;

end

save('statistics/significanceVAF.mat', 'pArray')