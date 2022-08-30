% sigDiffPaired.m
%  
% Calculate if the mean of population 1 is significantly larger than the
% mean of population 2. Critical values are those for 10 degrees of freedom
% (11 subjects).
function [T, sig5, sig25, sig1] = sigDiffPaired(pop1, pop2)
   
    T = mean(pop1-pop2)/(std(pop1-pop2)/sqrt(length(pop1)));

    if T > 2.764
        sig1 = true;
        sig25 = true;
        sig5 = true;
        
    elseif T > 2.228
        sig1 = false;
        sig25 = true;
        sig5 = true;
    
    elseif T > 1.812
        sig1 = false;
        sig25 = false;
        sig5 = true;
        
    else
        sig1 = false;
        sig25 = false;
        sig5 = false;
    end
    
end