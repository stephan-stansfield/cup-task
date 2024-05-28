% statistics_ANOVA
%
% Conduct Analysis of Variance on best-fit simulation data

load("best fit simulation/VAFarrays.mat");      % load best-fit VAF data

medianVAF = medianVAF(2:end,:);         % remove original input shaping results for ANOVA
length = numel(medianVAF);
y = reshape(medianVAF, [1, length]);    % flatten array (response variable is median VAF)

g_subj = horzcat(ones(1,10), 2*ones(1,10), 3*ones(1,10), 4*ones(1,10),...
    5*ones(1,10), 6*ones(1,10), 7*ones(1,10), 8*ones(1,10), 9*ones(1,10),...
    10*ones(1,10), 11*ones(1,10));

sim_vec = ["MM", "SM", "FM", "RB", "NI"];
g_sim = repmat(sim_vec,[1,22]);
% g_sim = horzcat(sim_array,sim_array,sim_array,sim_array,sim_array,...
%     sim_array,sim_array,sim_array,sim_array,sim_array,sim_array);

FF_vec = ["NF", "NF", "NF", "NF", "NF", "FF", "FF", "FF", "FF", "FF"];
g_FF = repmat(FF_vec,[1,11]);

[p,tbl,stats] = anovan(y,{g_subj,g_sim,g_FF})

%%
% Conduct ANOVA just on subset of simulations that include feedforward

y_FF = y(56:end);

g_subj_FF = horzcat(ones(1,5), 2*ones(1,5), 3*ones(1,5), 4*ones(1,5),...
    5*ones(1,5), 6*ones(1,5), 7*ones(1,5), 8*ones(1,5), 9*ones(1,5),...
    10*ones(1,5), 11*ones(1,5));

g_sim_FF = g_sim(1:55);

[p_FF,tbl_FF,stats_FF] = anovan(y_FF,{g_subj_FF,g_sim_FF})
