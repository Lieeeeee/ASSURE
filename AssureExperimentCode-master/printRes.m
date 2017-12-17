function [ meanVals ] = printRes( res,imgCell )
%PRINTRES Summary of this function goes here
%   Detailed explanation goes here

meanVals = [];
gmax = ['dice(MAX,est(MAX)) ' sprintf('%.3f ', res.max_dice)];
gmin = ['dice(MIN,est(MIN)) ' sprintf('%.3f ', res.min_dice)];
gvar = ['dice(VAR,est(VAR)) ' sprintf('%.3f ', res.var_dice)];
meanVals = [meanVals,mean(res.max_dice),mean(res.min_dice),mean(res.var_dice)];

%g_ubound = ['min_area_bound ' sprintf('%.3f ', res.var_percentage_range_gt(2,:))];
%g_ubounde_stimated = ['min_area_bound_estimated ' sprintf('%.3f ', res.var_percentage_range(2,:))];
g_ubound = ['possible ' sprintf('%.3f ', res.var_range_gt(2,:))];
g_ubounde_stimated = ['possible_estimated ' sprintf('%.3f ', res.var_range(2,:))];
g_ubound_diff = ['percent_change ' sprintf('%.3f ', (res.var_range_gt(2,:)-res.var_range(2,:))./res.var_range_gt(2,:))];
meanVals = [meanVals,mean((res.var_range_gt(2,:)-res.var_range(2,:))./res.var_range_gt(2,:))];

%g_lbound = ['max_area_bound ' sprintf('%.3f ', res.var_percentage_range_gt(1,:))];
%g_lbounde_stimated = ['max_area_bound_estimated ' sprintf('%.3f ', res.var_percentage_range(1,:))];
g_lbound = ['consensus ' sprintf('%.3f ', res.var_range_gt(1,:))];
g_lbounde_stimated = ['consensus_estimated ' sprintf('%.3f ', res.var_range(1,:))];
g_lbound_diff = ['percent_change ' sprintf('%.3f ', (res.var_range_gt(1,:)-res.var_range(1,:)) ./ res.var_range_gt(1,:))];
meanVals = [meanVals,mean((res.var_range_gt(1,:)-res.var_range(1,:))./res.var_range_gt(1,:))];

g_var_vol = ['var_volume ' sprintf('%.3f ', res.var_volGT)];
g_var_vol_estimated = ['var_volume_estimate ' sprintf('%.3f ', res.var_vol)];
g_var_vol_diff = ['percent_change ' sprintf('%.3f ', (res.var_volGT-res.var_vol) ./ res.var_volGT)];
meanVals = [meanVals,mean((res.var_volGT-res.var_vol) ./ res.var_volGT)];

%r_ubound = ['min_recist_bound ' sprintf('%.3f ', res.recist_percentage_range_gt(1,:))];
%r_ubounde_stimated = ['min_recist_bound_estimated ' sprintf('%.3f ', res.recist_percentage_range(1,:))];
r_ubound = ['recist_consensus ' sprintf('%.3f ', res.recist_range_gt(1,:))];
r_ubounde_stimated = ['recist_consensus_estimated ' sprintf('%.3f ', res.recist_range(1,:))];
r_ubound_diff = ['percent_change ' sprintf('%.3f ', (res.recist_range_gt(1,:)-res.recist_range(1,:))./res.recist_range_gt(1,:))];
meanVals = [meanVals,mean((res.recist_range_gt(1,:)-res.recist_range(1,:))./res.recist_range_gt(1,:))];

%r_lbound = ['max_recist_bound ' sprintf('%.3f ', res.recist_percentage_range_gt(2,:))];
%r_lbounde_stimated = ['max_recist_bound_estimated ' sprintf('%.3f ', res.recist_percentage_range(2,:))];
r_lbound = ['recist_possible ' sprintf('%.3f ', res.recist_range_gt(2,:))];
r_lbounde_stimated = ['recist_possible_estimated ' sprintf('%.3f ', res.recist_range(2,:))];
r_lbound_diff = ['percent_change  ' sprintf('%.3f ', (res.recist_range_gt(2,:)-res.recist_range(2,:))./res.recist_range_gt(2,:))];
meanVals = [meanVals,mean((res.recist_range_gt(2,:)-res.recist_range(2,:))./res.recist_range_gt(2,:))];

g_vol_div_area = ['var-vol-area/mean-area(N) ' sprintf('%.3f ', res.var_volGT./res.meanPixAreas')];

g_l_percent = ['vol_percent_low_diff ' sprintf('%.3f ', (res.var_percentage_range_gt(1,:)))];
g_u_percent = ['vol_percent_upper_diff ' sprintf('%.3f ', (res.var_percentage_range_gt(2,:)))];
g_l_percent_estimated = ['vol_percent_low_diff_estimated ' sprintf('%.3f ', (res.var_percentage_range(1,:)))];
g_u_percent_estimated = ['vol_percent_upper_diff_estimated ' sprintf('%.3f ', (res.var_percentage_range(2,:)))];


r_l_percent = ['recist_percent_low_diff ' sprintf('%.3f ', (res.recist_percentage_range_gt(1,:)))];
r_u_percent = ['recist_percent_upper_diff ' sprintf('%.3f ', (res.recist_percentage_range_gt(2,:)))];
r_l_percent_estimated = ['recist_percent_low_diff_estimated ' sprintf('%.3f ', (res.recist_percentage_range(1,:)))];
r_u_percent_estimated = ['recist_percent_upper_diff_estimated ' sprintf('%.3f ', (res.recist_percentage_range(2,:)))];



fprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n\n', gmax,gmin,gvar,g_ubound,g_ubounde_stimated,g_ubound_diff,g_lbound,g_lbounde_stimated,g_lbound_diff, ...
    r_ubound,r_ubounde_stimated,r_ubound_diff,r_lbound,r_lbounde_stimated,r_lbound_diff,g_l_percent,g_u_percent,g_l_percent_estimated,g_u_percent_estimated, ...
    r_l_percent,r_u_percent,r_l_percent_estimated,r_u_percent_estimated,g_var_vol,g_var_vol_estimated,g_var_vol_diff);


end

