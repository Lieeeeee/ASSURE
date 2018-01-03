% run a genetic algorithm on the parameter space of snakes Options
% the parameters we are optimizing on:
%   iterations
%   Wline
%   Wedge
%   Wterm
%   alpha
%   beta
%   delta
%   kappa

nvars = 8;
% naiive try - without constraining the problem and without setting the
% ranges of the params
% [x, fval] = ga(@getPveBasedLoss, 2*nvars);

% setting the options
% initial population
init_it_in = 2; init_wline_in = 0.04; init_wedge_in = 2.0; init_wterm_in = 0.01; init_alpha_in = 0.2;
init_beta_in = 0.2; init_delta_in = -0.5; init_kappa_in = 2; 
init_it_out = 3; init_wline_out = -0.04; init_wedge_out = -2.0; init_wterm_out = -0.01; init_alpha_out = 0.2;
init_beta_out = 0.2; init_delta_out = 0.5; init_kappa_out = 2; 
initial_population = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
    init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out];
iterations = 50*nvars;
% creating options object
options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'InitialPopulationMatrix', initial_population, 'Display', 'iter', 'MaxGenerations', iterations);

% calling ga still without setting constraints
[x, fval] = ga(@getPveBasedLoss, 2*nvars, [],[],[],[],[],[],[], options);