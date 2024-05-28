function sigdiff(x1,x2,s1,s2,n1,n2)
    
    % Critical values (two-tailed, large sample size)
    C1p = 2.576;
    C2p = 2.326;
    C5p = 1.96;

    Z = abs((x1 - x2)/sqrt((s1^2)/n1 + (s2^2)/n2));
    if Z > C1p
        disp('1%')
    elseif Z > C2p
        disp('2.5%')
    elseif Z > C5p
        disp('5%')
    else
        disp('Not significant')
    end
    
end