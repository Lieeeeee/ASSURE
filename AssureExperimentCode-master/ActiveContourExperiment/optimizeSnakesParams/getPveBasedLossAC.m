function [ loss ] = getPveBasedLossAC( params )
%GETPVEBASEDLOSSAC Summary of this function goes here
%   Detailed explanation goes here

%GETPVEBASEDLOSS this will be used as an objective function for
% parameter optimization. 
% holds an experiment, and then calculates and returns the percentage of error in the
% variability estimation. 
%   INPUT:
%       params - the parameters we are now optimizing on. as a vector of
%       length n_params
%   OUTPUT:
%       loss - the loss from the resulting experiment

    % load the images
    load('chapter3Env.mat');
    imgCell = brainImgCell3D;
    
    % prepare the parameters for the experiment
    % % lung
    % kernelSize = 2; Tlength = 3;
    % liver
    kernelSize = 4; Tlength = 3;
    SNAKE_RUN = 1;
    TRIVIAL_RUN = 0;
    SEGNUM_TO_USE = 0;

    % extract optimized parameters
    iterIn = round(params(1)); cIn = params(2); sIn = max(0, params(3)); 
    iterOut = round(params(4)); cOut = params(5); sOut = max(0, params(6)); 
    
    % create an Options struct for in and out
    OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
    OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
    
    % run experiment
    [res, ~] =  VariabilityExperiment.holdExperiment(imgCell, kernelSize, [], [], Tlength, TRIVIAL_RUN, SNAKE_RUN, ...
        SEGNUM_TO_USE, OptionsIn, OptionsOut);
    
    % compute loss and return
%     loss = mean(abs((res.var_volGT-res.var_vol) ./ res.var_volGT)) + mean(abs((res.var_range_gt(2,:)-res.var_range(2,:))./res.var_range_gt(2,:))) + ...
%         mean(abs((res.var_range_gt(1,:)-res.var_range(1,:)) ./ res.var_range_gt(1,:))) + (1 - min(res.var_dice));
    loss = - min(res.var_dice);

end
