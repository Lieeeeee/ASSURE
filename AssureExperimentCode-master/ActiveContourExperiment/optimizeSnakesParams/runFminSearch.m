%% use fminsearch in order to find a local minima of the Pve

% % default
% wline_default = 0.04; wedge_default = 2.0; wterm_default = 0.01;
% alpha_default = 0.2; beta_default = 0.2; delta_default = 0.1; kappa_default = 2;
% 
% % % lung
% % init_it_in = 2; init_wline_in = 0.04; init_wedge_in = 2.0; init_wterm_in = 0.01; init_alpha_in = 0.2;
% % init_beta_in = 0.2; init_delta_in = -0.5; init_kappa_in = 2; 
% % init_it_out = 3; init_wline_out = -0.04; init_wedge_out = -2.0; init_wterm_out = -0.01; init_alpha_out = 0.2;
% % init_beta_out = 0.2; init_delta_out = 0.5; init_kappa_out = 2; 
% % initial_population = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
% %     init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out];
% 
% % liver
% init_it_in = 2; init_wline_in = -0.5; init_wedge_in = wedge_default; init_wterm_in = wterm_default; init_alpha_in = alpha_default;
% init_beta_in = beta_default; init_delta_in = -0.05; init_kappa_in = kappa_default; 
% init_it_out = 3; init_wline_out = wline_default; init_wedge_out = wedge_default; init_wterm_out = 0.1; init_alpha_out = alpha_default;
% init_beta_out = beta_default; init_delta_out = 1; init_kappa_out = -2; 
% initial_population = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
%     init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out];
% 
% options = optimset('Display','iter');
% [x, fval] = fminsearch(@getPveBasedLoss, initial_population, options);
% save('AssureExperimentCode-master/ActiveContourExperiment/optimizeSnakesParams/fminSearchDice_liver_maxpenalty_without_dice.mat')

%% activeContours MATLAB implementation
% lung
initial_population = [4, -0.5, 1, 4, 0.1, 1]; 
% initial_population = [3, -4.4412, 1.3351, 4, 9.0265, 6.1125];
% options = optimset('Display','iter');
% [x, fval] = fminsearch(@getPveBasedLossAC, initial_population, options);
lb = [1, -Inf, 0, 1, -Inf, 0]; ub = [Inf, Inf, Inf, Inf, Inf, Inf];
% [x, fval] = fmincon(@getPveBasedLossAC,initial_population,[],[],[],[],lb,ub,[],options);

% GA
nvars = 6;
generations = 5*nvars;
% creating options object
options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', initial_population, 'Display', 'iter', 'MaxGenerations', generations);

% calling ga still without setting constraints
[x, fval] = ga(@getPveBasedLossAC, nvars, [],[],[],[], lb, ub, [], options);
