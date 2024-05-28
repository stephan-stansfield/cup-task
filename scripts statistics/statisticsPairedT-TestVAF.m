% statistics_paired_ttest
%
% Conduct dependent t-test for paired samples. See: 
% https://en.wikipedia.org/wiki/Student%27s_t-test#Dependent_t-test_for_paired_samples

load("best fit simulation/VAFarrays.mat");      % load best-fit VAF data

n = 11; % number of subjects
nSim = 5;
difference = zeros(nSim,n); % difference between subject mean VAFs
X = zeros(1,nSim); % mean of differences between subject pairs
s = zeros(1,nSim); % std dev of differences between subject pairs
% u0 = mean(medianVAF(1,:)); % mean of VAF across subjects for original input shaping
u0 = 0;
sigArray = strings(1,nSim);

for sim = 1:nSim
    difference(sim,:) = meanVAF(sim+1,:) - meanVAF(1,:); % compare each simulation type to nominal input shaping
    X(sim) = mean(difference(sim,:));
    s(sim) = std(difference(sim,:));
    sigArray(sim) = sigdiff(X(sim),u0,s(sim),n);
end


function significance = sigdiff(x,u,s,n)
    % Critical values (one-tailed, dof = 10)
    C1p = 2.764;
    C25p = 2.228;
    C5p = 1.812;

    t = (x-u)/(s/sqrt(n))
    if t > C1p
        significance = '**';
    elseif t > C25p
        significance = '*/';
    elseif t > C5p
        significance = '*';
    else
        significance = '-';
    end

end
