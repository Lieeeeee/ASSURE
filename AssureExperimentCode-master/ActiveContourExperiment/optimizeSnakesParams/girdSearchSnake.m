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
%     sigma1_range = [3:15];
    Wline_range = [-2:0.2:2];
    Wedge_range = [-2:0.2:2];
    Wterm_range = [-2:0.2:2];
%     sigma2_range = [3:15];
%     mu_range = [0:0.01:0.5];
%     GIter_range = [1:10];
%     sigma3_range = [3:15];
    alpha_range = [-2:0.2:2];
    beta_range = [-2:0.2:2];
    delta_range = [-2:0.2:2];
    kappa_range = [-2:0.2:2];

    % get grid
    [iterationsIn, WlineIn, WedgeIn, WtermIn, alphaIn, betaIn, deltaIn, kappaIn, ... 
        iterationsOut, WlineOut, WedgeOut, WtermOut, alphaOut, betaOut, deltaOut, kappaOut] = ...
        ndgrid(iter_range, Wline_range, Wedge_range, Wterm_range, alpha_range, beta_range, delta_range, kappa_range, ...
        iter_range, Wline_range, Wedge_range, Wterm_range, alpha_range, beta_range, delta_range, kappa_range);
%     even this is too big
%     [iterationsIn, WlineIn, WedgeIn, WtermIn, alphaIn, betaIn, deltaIn, kappaIn] = ...
%         ndgrid(iter_range, Wline_range, Wedge_range, Wterm_range, alpha_range, beta_range, delta_range, kappa_range);

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
    [~] = arrayfun(@(iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn, ...
        iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut, ...
        imgCell, experiment_folder, outputFile, kernelSize, Tlength) runSingleGridSearchIteration ...
        (iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn, ...
        iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut, imgCell, experiment_folder, outputFile, kernelSize, Tlength), ...
        iterationsIn, WlineIn, WedgeIn, WtermIn, alphaIn, betaIn, deltaIn, kappaIn, ... 
        iterationsOut, WlineOut, WedgeOut, WtermOut, alphaOut, betaOut, deltaOut, kappaOut, ...
        imgCell, experiment_folder, outputFileAll, kernelSize, Tlength);
   
end

