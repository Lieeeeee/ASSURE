classdef ActiveContourOptions
    %ACTIVECONTOUROPTIONS Summary of this class goes here 
    %     n_iterations - number of iterations
    %     contraction_param_out, contraction_param_in - Tendency of the
    %       contour to grow outwards or shrink inwards. Positive values bias
    %       the contour to shrink inwards (contract).
    %     smooth_param_out, smooth_param_in - Degree of smoothness or
    %       regularity of the boundaries of the segmented regions. Higher
    %       values produce smoother region boundaries.
    
    properties
    end
    
    methods(Static)
        function [Options] = getEmptyOptions()
            Options = struct;
        end
        
        function [Options] = getSpecifiedOptions(iter, contration, smoothness)
            Options = ActiveContourOptions.getEmptyOptions();
            Options.Iterations = iter; 
            Options.Contraction = contration;
            Options.Smoothness = smoothness;
        end
        
        function [OptionsIn, OptionsOut] = getLungOptions()
%             params = [3, -4.4412, 1.3351, 4, 9.0265, 6.1125]; % GA result
            params = [4, 1.2031, 0.8832, 4, -0.3484, 0.1173]; % fminSearch result
%             params = [4, 1, 1, 4, -0.5, 0.1];
            iterIn = params(1); cIn = params(2); sIn = params(3); 
            OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
            
            iterOut = params(4); cOut = params(5); sOut = params(6); 
            OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
        end
        
        function [OptionsIn, OptionsOut] = getLiverOptions()
            params = [4, 1.1364, 1.1158, 4, -0.2592, 0.1026]; % fminSearch result
%             params = [4, 1, 1, 4, -0.5, 0.1];
            iterIn = params(1); cIn = params(2); sIn = params(3); 
            OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
            
            iterOut = params(4); cOut = params(5); sOut = params(6); 
            OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
        end
        
        function [OptionsIn, OptionsOut] = getKidneyOptions()
            params = [4, 0.8629, 1.3791, 4, -0.2840, 0.1197]; % fminSearch result
%             params = [4, 1, 1, 4, -0.5, 0.1];
            iterIn = params(1); cIn = params(2); sIn = params(3); 
            OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
            
            iterOut = params(4); cOut = params(5); sOut = params(6); 
            OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
        end
        
        function [OptionsIn, OptionsOut] = getBrainOptions()

            params = [4, 1.3343, 0.1312, 4, -0.7082, 0.0933]; % fminSearch result
%             params = [4, 1, 1, 4, -0.5, 0.1];
            iterIn = params(1); cIn = params(2); sIn = params(3); 
            OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
            
            iterOut = params(4); cOut = params(5); sOut = params(6); 
            OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
        end
%         
%         function [OptionsStruct] = getAllOptions(DEBUG)
%             OptionsStruct = struct;
%             
%             % lung
%             [OptionsIn, OptionsOut] = SnakeOptions.getLungOptions(DEBUG);
%             OptionsStruct.LungOptions = {OptionsIn, OptionsOut};
%             % liver
%             [OptionsIn, OptionsOut] = SnakeOptions.getLiverOptions(DEBUG);
%             OptionsStruct.LiverOptions = {OptionsIn, OptionsOut};
%             % kidney
%             [OptionsIn, OptionsOut] = SnakeOptions.getKidneyOptions(DEBUG);
%             OptionsStruct.KidneyOptions = {OptionsIn, OptionsOut};
%             % brain
%             [OptionsIn, OptionsOut] = SnakeOptions.getBrainOptions(DEBUG);
%             OptionsStruct.BrainOptions = {OptionsIn, OptionsOut};
%         end
    
    end
end

