% statistics_paired_ttest_peaks
%
% Conduct dependent 2-tailed t-test for paired samples, applied to human
% data vs. simulated peak ratios. Test if simulated median peak ratios 
% paired by subject siginificantly differ from human data. See: 
% https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples

load("peakRatioArrays.mat");      % load best-fit peak ratio data

n = 11;                     % number of subjects
nSim = 6;                   % number of simulations
difference = zeros(nSim,n); % difference between subject mean peak ratios
X = zeros(1,nSim);          % mean of differences between subject pairs
s = zeros(1,nSim);          % std dev of differences between subject pairs
u0 = 0;
sigArray = strings(nSim,1);
tArray = zeros(nSim,1);

for sim = 1:nSim
    difference(sim,:) = meanBySubject(sim+1,:) - meanBySubject(1,:); % compare each simulation type to human data
    X(sim) = mean(difference(sim,:));
    s(sim) = std(difference(sim,:));
    [sigArray(sim), tArray(sim)] = sigdiff(X(sim),u0,s(sim),n);
end


function [significance, t] = sigdiff(x,u,s,n)
    % Critical values (two-tailed, dof = 10)
    C1p = 3.169;
    C2p = 2.764;
    C5p = 2.228;

    t = abs(x-u)/(s/sqrt(n))
    if t > C1p
        significance = '**';
    elseif t > C2p
        significance = '*/';
    elseif t > C5p
        significance = '*';
    else
        significance = '-';
    end

end
