% use fminsearch in order to find a local minima of the Pve

% previous good params, resulted in loss of 0.0523!
initial_population = [2.5451    0.0334    2.1619    0.0088    0.1600    0.1711   -0.5740    2.3761    3.5504   -0.0586   -2.0904   -0.0071    0.1349    0.2608    0.4343    3.0997];

options = optimset('Display','iter');

[x, fval] = fminsearch(@getPveBasedLoss, initial_population, options);