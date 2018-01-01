function [] = girdSearchSnake( imgCell, experiment_folder, kernelSize, Tlength )
%GIRDSEARCHSNAKE performs a grid search on the parameter space of the
%snakes active contour function in order to find the best parameters for
%the given structure
%   INPUT:
%       imgCell - the image cell to work on in these experiments
%       experiment_folder - the folder to write the results to
%       outputFileMat - the file to write the results to
%       kernelSize, Tlength - parameters for the 'holdExperiment' function

    % define ranges
    iter_range = [1:10];
    sigma1_range = [3:15];
    Wline_range = [-5:0.01:5];
    Wedge_range = [-5:0.01:5];
    Wterm_range = [-5:0.01:5];
    sigma2_range = [3:15];
    mu_range = [0:0.01:0.5];
    GIter_range = [1:10];
    sigma3_range = [3:15];
    alpha_range = [-5:0.01:5];
    beta_range = [-5:0.01:5];
    delta_range = [-5:0.01:5];
    kappa_range = [-5:0.01:5];

    % get grid
    [iterationsIn, sigma1In, WlineIn, WedgeIn, WtermIn, sigma2In, muIn, GIterationsIn, sigma3In, alphaIn, betaIn, deltaIn, kappaIn, ... 
        iterationsOut, sigma1Out, WlineOut, WedgeOut, WtermOut, sigma2Out, muOut, GIterationsOut, sigma3Out, alphaOut, betaOut, deltaOut, kappaOut] = ...
        ndgrid(iter_range, sigma1_range, Wline_range, Wedge_range, Wterm_range, sigma2_range, mu_range, GIter_range, sigma3_range, alpha_range, ...
        beta_range, delta_range, kappa_range, ...
        iter_range, sigma1_range, Wline_range, Wedge_range, Wterm_range, sigma2_range, mu_range, GIter_range, sigma3_range, alpha_range, beta_range, ...
        delta_range, kappa_range);
    
    % open file for all the experiments to write to
    outputFileAll = [experiment_folder '/gridSearch.csv'];
    fid = fopen(outputFileAll, 'w');
    col_defs = {'VerboseIn','nPointsIn','WlineIn','WedgeIn','WtermIn','Sigma1In','Sigma2In',...
        'AlphaIn','BetaIn','DeltaIn','GammaIn','KappaIn','IterationsIn','GIterationsIn','MuIn','Sigma3In', ...
        'VerboseOut','nPointsOut','WlineOut','WedgeOut','WtermOut','Sigma1Out','Sigma2Out',...
        'AlphaOut','BetaOut','DeltaOut','GammaOut','KappaOut','IterationsOut','GIterationsOut','MuOut','Sigma3Out', ...
        'possible_dice', 'concensus_dice', 'var_dice', 'possible_gt', 'possible_est', 'possible_per_change', ...
        'consensus_gt', 'consensus_est', 'consensus_per_change', 'var_gt', 'var_est', 'var_per_change'};
    fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', col_defs{1, 1:end});
    fclose(fid);
    
    % run the snakes active contour on all the parameters in the grid
    [~] = arrayfun(@(iterIn, sig1In, WlIn, WeIn, WtIn, sig2In, mIn, GIterIn, sig3In, alIn, beIn, delIn, kapIn, ...
        iterOut, sig1Out, WlOut, WeOut, WtOut, sig2Out, mOut, GIterOut, sig3Out, alOut, beOut, delOut, kapOut, ...
        imgCell, experiment_folder, outputFile, kernelSize, Tlength) runSingleGridSearchIteration ...
        (iterIn, sig1In, WlIn, WeIn, WtIn, sig2In, mIn, GIterIn, sig3In, alIn, beIn, delIn, kapIn, ...
        iterOut, sig1Out, WlOut, WeOut, WtOut, sig2Out, mOut, GIterOut, sig3Out, alOut, beOut, delOut, kapOut, imgCell, experiment_folder, outputFile, kernelSize, Tlength), ...
        iterationsIn, sigma1In, WlineIn, WedgeIn, WtermIn, sigma2In, muIn, GIterationsIn, sigma3In, alphaIn, betaIn, deltaIn, kappaIn, ... 
        iterationsOut, sigma1Out, WlineOut, WedgeOut, WtermOut, sigma2Out, muOut, GIterationsOut, sigma3Out, alphaOut, betaOut, deltaOut, kappaOut, ...
        imgCell, experiment_folder, outputFileAll, kernelSize, Tlength);
   
end

