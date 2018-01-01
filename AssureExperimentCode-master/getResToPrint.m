function [ resToWrite ] = getResToPrint( res )
%GETRESTOPRINT prepare a mat to later write to a file
%   INPUT:
%       res - the results struct

    possible_dice = res.max_dice; 
    concensus_dice = res.min_dice;
    var_dice = res.var_dice;
    possible_gt = res.var_range_gt(2,:);
    possible_est = res.var_range(2,:);
    possible_per_change = (res.var_range_gt(2,:)-res.var_range(2,:))./res.var_range_gt(2,:);
    consensus_gt = res.var_range_gt(1,:);
    consensus_est = res.var_range(1,:);
    consensus_per_change = (res.var_range_gt(1,:)-res.var_range(1,:)) ./ res.var_range_gt(1,:);
    var_gt = res.var_volGT;
    var_est = res.var_vol;
    var_per_change = (res.var_volGT-res.var_vol) ./ res.var_volGT;

%     recist_consensus_gt = res.recist_range_gt(1,:);
%     recist_consensus_est = res.recist_range(1,:);
%     recist_consensus_diff = (res.recist_range_gt(1,:)-res.recist_range(1,:))./res.recist_range_gt(1,:);
%     recist_possible_gt = res.recist_range_gt(2,:);
%     recist_possible_est = res.recist_range(2,:);
%     recist_possible_diff = (res.recist_range_gt(2,:)-res.recist_range(2,:))./res.recist_range_gt(2,:);
    
    resToWrite = [possible_dice, concensus_dice, var_dice, ...
        possible_gt, possible_est, possible_per_change, ...
        consensus_gt, consensus_est, consensus_per_change, ...
        var_gt, var_est, var_per_change];
end

