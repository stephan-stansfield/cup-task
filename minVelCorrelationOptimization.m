% minVelCorrelationOptimization.m

Bmin = 6;
Bmax = 10;
Kmin = 75;
Kmax = 125;

function [algorithm, lb, ub, maxeval] = NLopt()
    
    opt.algorithm = NLOPT_GN_CRS2_LM;
    opt.lower_bounds = [Bmin, Kmin];
    opt.upper_bounds = [Bmax, Kmax];
    init_guess = [Bmin + (Bmax - Bmin)/2, Kmin + (Kmax - Kmin)/2;
        
    opt.max_objective = @(x) corrObjFunc(x);
    
    [xopt, fmin, ~] = nlopt_optimize(opt, init_guess);
    
end

