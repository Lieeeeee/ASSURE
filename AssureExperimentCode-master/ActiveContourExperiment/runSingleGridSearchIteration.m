function [] = runSingleGridSearchIteration( iterIn, sig1In, WlIn, WeIn, WtIn, sig2In, mIn, GIterIn, sig3In, alIn, beIn, delIn, kapIn, ...
        iterOut, sig1Out, WlOut, WeOut, WtOut, sig2Out, mOut, GIterOut, sig3Out, alOut, beOut, delOut, kapOut, ...
        imgCell, experiment_folder, outputFile, kernelSize, Tlength )
%runSingleGridSearchIteration creates a snakeOptions object with the given
%parameters and runs an experiment with them. 
%   INPUT: 
%       gets all the parameters of Options we want to test twice: for in
%       and for out
%       imgCell - the image cell to work on in these experiments
%       experiment_folder - the folder to write the results to -- don't
%       want to write nothing so maybe unneccessary
%       outputFile - the file to write the results to
%       kernelSize, Tlength - parameters for the 'holdExperiment' function
    
    SNAKE_RUN = 2;
    TRIVIAL_RUN = 0;
    SEGNUM_TO_USE = 0;

    % create an Options struct for in and out
    OptionsIn = SnakeOptions.getSpecifiedOptions(false, iterIn, sig1In, WlIn, WeIn, WtIn, sig2In, mIn, GIterIn, sig3In, alIn, beIn, delIn, kapIn);
    OptionsOut = SnakeOptions.getSpecifiedOptions(false, iterOut, sig1Out, WlOut, WeOut, WtOut, sig2Out, mOut, GIterOut, sig3Out, alOut, beOut, delOut, kapOut);
    
    % run experiment
    [res, ~] =  VariabilityExperiment.holdExperiment(imgCell, kernelSize, [], [], Tlength, TRIVIAL_RUN, SNAKE_RUN, ...
        SEGNUM_TO_USE, OptionsIn, OptionsOut);
    
    % prepare ouput to write
    optionsToWrite = [cell2mat(struct2cell(OptionsInLung)); cell2mat(struct2cell(OptionsOutLung))]';
    resToWrite = getResToPrint(res);
    toWrite = [optionsToWrite, resToWrite];
        
    % write results to file
    prev = csvread(outputFile);
    csvwrite(outputFile, [prev; toWrite]);
end

