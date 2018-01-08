%% use fminsearch in order to find a local minima of the Pve

% previous good params
% initial_population = [2.9388    0.0377    1.5447    0.0108    0.1367    0.1681   -0.2971    2.6419    3.0931   -0.0576   -2.6526   -0.0094    0.1291    0.1801    0.5126    3.3646];
init_it_in = 2; init_wline_in = 0.04; init_wedge_in = 2.0; init_wterm_in = 0.01; init_alpha_in = 0.2;
init_beta_in = 0.2; init_delta_in = -0.5; init_kappa_in = 2; 
init_it_out = 3; init_wline_out = -0.04; init_wedge_out = -2.0; init_wterm_out = -0.01; init_alpha_out = 0.2;
init_beta_out = 0.2; init_delta_out = 0.5; init_kappa_out = 2; 
initial_population = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
    init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out];

options = optimset('Display','iter');

[x, fval] = fminsearch(@getPveBasedLoss, initial_population, options);
save('AssureExperimentCode-master/ActiveContourExperiment/optimizeSnakesParams/fminSearchDice.mat')