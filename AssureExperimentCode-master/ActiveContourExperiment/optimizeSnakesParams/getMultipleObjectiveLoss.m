function [ loss ] = getMultipleObjectiveLoss( params )
%GETMULTIPLEOBJECTIVELOSS this will be used as an objective function for
% parameter optimization when we want to achieve multiple objectives:
%   minimum percent change in all consensus, possible, and variability;
%   maximum dice coefficient for the variability area
% holds an experiment, and then calculates and returns the percentage of
% error.
%   INPUT:
%       params - the parameters we are now optimizing on. as a vector of
%       length n_params
%   OUTPUT:
%       loss - the loss from the resulting experiment

    % load the images
    load('chapter3Env.mat');
    imgCell = lungImgCell3D;
    
    % prepare the parameters for the experiment
    kernelSize = 2; Tlength = 3;
    SNAKE_RUN = 2;
    TRIVIAL_RUN = 0;
    SEGNUM_TO_USE = 0;
    DEBUG = false;

    % extract optimized parameters
    iterIn = params(1); WlIn = params(2); WeIn = params(3); WtIn = params(4); alIn = params(5); beIn = params(6); delIn = params(7); kapIn = params(8);
    iterOut = params(9); WlOut = params(10); WeOut = params(11); WtOut = params(12); alOut = params(13); beOut = params(14); delOut = params(15); kapOut = params(16);
    
    % create an Options struct for in and out
    OptionsIn = SnakeOptions.getSpecifiedOptions(DEBUG, iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn);
    OptionsOut = SnakeOptions.getSpecifiedOptions(DEBUG, iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut);
    
    % run experiment
    [res, ~] =  VariabilityExperiment.holdExperiment(imgCell, kernelSize, [], [], Tlength, TRIVIAL_RUN, SNAKE_RUN, ...
        SEGNUM_TO_USE, OptionsIn, OptionsOut);
    
    % compute loss and return
    loss = zeros(4,1);
    loss(1) = mean(abs((res.var_range_gt(2,:)-res.var_range(2,:))./res.var_range_gt(2,:)));
    loss(2) = mean(abs((res.var_range_gt(1,:)-res.var_range(1,:)) ./ res.var_range_gt(1,:)));
    loss(3) = mean(abs((res.var_volGT-res.var_vol) ./ res.var_volGT));
    loss(4) = -mean(res.var_dice);

end

