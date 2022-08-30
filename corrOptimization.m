% corrOptimization.m

Bmin = 6;
Bmax = 10;
Kmin = 75;
Kmax = 125;

plotResults = true;

opt.maxeval = 1000;
opt.xtol_rel = 0.00001;
opt.verbose = 1;
% opt.algorithm = NLOPT_GN_CRS2_LM;
opt.algorithm = NLOPT_GN_DIRECT_L;
opt.lower_bounds = [Bmin, Kmin];
opt.upper_bounds = [Bmax, Kmax];
init_guess = [Bmin + (Bmax - Bmin)/2, Kmin + (Kmax - Kmin)/2];

opt.max_objective = @(x) corrObjFunc(x,false);

[xopt, fmin, ~] = nlopt_optimize(opt, init_guess);

disp('Best fit values:')
optB = xopt(1)
optK = xopt(2)

if plotResults
    
   corrObjFunc(xopt, plotResults)
    
end

