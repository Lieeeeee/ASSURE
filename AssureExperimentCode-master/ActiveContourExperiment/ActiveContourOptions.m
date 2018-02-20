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
%             params = [3, -4.4412, 1.3351, 4, 9.0265, 6.1125];
            params = [4, 1, 1, 4, -0.5, 0.1];
            iterIn = params(1); cIn = params(2); sIn = params(3); 
            OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
            
            iterOut = params(4); cOut = params(5); sOut = params(6); 
            OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);
        end
        
%         function [OptionsIn, OptionsOut] = getLiverOptions(DEBUG)
%             params = [1.6645   -0.6115    3.8542    0.0160    0.3467    0.0264   -0.1131    1.0516    3.8982    0.0516   -6.8314    0.3190   -0.0825    0.4198    0.6603   -2.2152]; % multiple objectives with 1-dice
% %             params = [1.4276   -0.3007    2.4772    0.0108    0.2288    0.2311   -0.0613    1.7769    3.5745    0.0341    1.8300    0.1016    0.2066    0.2252    1.0686   -1.3575]; % multiple objectives with 1/dice and max
% %             params = [1.4527   -0.4702    1.9121    0.0076    0.2302    0.2271   -0.0794    2.5498    3.9983    0.0258    2.7338    0.0626    0.1784    0.1640    0.8441   -1.2624];  % max penalty
% %             params = [1.8582   -0.5917    1.7512    0.0101    0.2161    0.2068   -0.0473    2.0106    3.3873    0.0437    2.4288    0.0999    0.1396    0.1691    0.9379   -2.2118]; % max penalty no dice
% %             init_it_in = 2; init_wline_in = -0.5; init_wedge_in = wedge_default; init_wterm_in = wterm_default; init_alpha_in = alpha_default;
% %             init_beta_in = beta_default; init_delta_in = -0.05; init_kappa_in = kappa_default; 
% %             init_it_out = 3; init_wline_out = wline_default; init_wedge_out = wedge_default; init_wterm_out = 0.1; init_alpha_out = alpha_default;
% %             init_beta_out = beta_default; init_delta_out = 1; init_kappa_out = -2; 
% %             params = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
% %                 init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out]; % original parameters
%             iterIn = params(1); WlIn = params(2); WeIn = params(3); WtIn = params(4); alIn = params(5); beIn = params(6); delIn = params(7); kapIn = params(8);
%             OptionsIn = SnakeOptions.getSpecifiedOptions(DEBUG, iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn);
%             
%             iterOut = params(9); WlOut = params(10); WeOut = params(11); WtOut = params(12); alOut = params(13); beOut = params(14); delOut = params(15); kapOut = params(16);
%             OptionsOut = SnakeOptions.getSpecifiedOptions(DEBUG, iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut);
%         end
%         
%         function [OptionsIn, OptionsOut] = getKidneyOptions(DEBUG)
%             OptionsIn = struct;
%             OptionsIn.Iterations = 2;
%             OptionsIn.Delta = -0.5;
% 
%             OptionsOut = struct;
%             OptionsOut.Iterations = 2;
%             OptionsOut.Delta = 1;
%             OptionsOut.Kappa = -2;
%             OptionsOut.Wterm = 0.1;
%             
%             if(DEBUG)
%                 OptionsIn.Verbose = true;
%                 OptionsOut.Verbose = true;
%             end
%         end
%         
%         function [OptionsIn, OptionsOut] = getBrainOptions(DEBUG)
%             OptionsIn = struct;
%             OptionsIn.Iterations = 5;
%             OptionsIn.Delta = -0.5;
% 
%             OptionsOut = struct;
%             OptionsOut.Iterations = 2;
%             OptionsOut.Delta = 1;
%             OptionsOut.Kappa = -2;
%             OptionsOut.Wterm = 0.1;
%             
%             if(DEBUG)
%                 OptionsIn.Verbose = true;
%                 OptionsOut.Verbose = true;
%             end
%         end
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

