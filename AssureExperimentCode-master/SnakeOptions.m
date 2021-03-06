classdef SnakeOptions
    %SNAKEOPTIONS class for defining the options for the snake algorithm
    % options (general),
    %  Option.Verbose : If true show important images, default false
    %  Options.nPoints : Number of contour points, default 100
    %  Options.Gamma : Time step, default 1
    %  Options.Iterations : Number of iterations, default 100
    %
    % options (Image Edge Energy / Image force))
    %  Options.Sigma1 : Sigma used to calculate image derivatives, default 10
    %  Options.Wline : Attraction to lines, if negative to black lines otherwise white
    %                    lines , default 0.04
    %  Options.Wedge : Attraction to edges, default 2.0
    %  Options.Wterm : Attraction to terminations of lines (end points) and
    %                    corners, default 0.01
    %  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
    %                    image (which gives the image force), default 20
    %
    % options (Gradient Vector Flow)
    %  Options.Mu : Trade of between real edge vectors, and noise vectors,
    %                default 0.2. (Warning setting this to high >0.5 gives
    %                an instable Vector Flow)
    %  Options.GIterations : Number of GVF iterations, default 0
    %  Options.Sigma3 : Sigma used to calculate the laplacian in GVF, default 1.0
    %
    % options (Snake)
    %  Options.Alpha : Membrame energy  (first order), default 0.2
    %  Options.Beta : Thin plate energy (second order), default 0.2
    %  Options.Delta : Baloon force, default 0.1
    %  Options.Kappa : Weight of external image force, default 2
    
    properties
    end
    
    methods(Static)
        function [Options] = getEmptyOptions(DEBUG)
            Options = struct;
            if(DEBUG)
                Options.Verbose = true;
            end
        end
        
        function [Options] = getSpecifiedOptions(DEBUG, iter, Wl, We, Wt, al, be, del, kap)
            Options = SnakeOptions.getEmptyOptions(DEBUG);
            Options.Iterations = iter; 
            Options.Wline = Wl;
            Options.Wedge = We;
            Options.Wterm = Wt;
            Options.Alpha = al;
            Options.Beta = be;
            Options.Delta = del;
            Options.Kappa = kap;
        end
        
        function [OptionsIn, OptionsOut] = getLungOptions(DEBUG)
%             params = [2.5451    0.0334    2.1619    0.0088    0.1600 0.1711   -0.5740    2.3761    3.5504   -0.0586   -2.0904 -0.0071    0.1349    0.2608    0.4343    3.0997]; % volume error
%             params = [2.9388    0.0377    1.5447    0.0108    0.1367    0.1681   -0.2971    2.6419    3.0931   -0.0576   -2.6526   -0.0094    0.1291    0.1801    0.5126    3.3646]; % multiple objectives
%             params = [2.7192    0.0539    1.4224    0.0305    0.0481   -0.0115   -0.7977    1.0008    3.9976   -0.0342    0.1349   -0.0114    0.2256    0.2467    0.8433    0.2744]; % dice objective
            params = [2.5289    0.0452    1.3843    0.0081    0.2009    0.0939   -0.4978    1.9144    3.4726   -0.0506   -2.3274   -0.0077    0.1705    0.1364    0.4721    3.6061]; % multiple objectives with 1/dice
            iterIn = params(1); WlIn = params(2); WeIn = params(3); WtIn = params(4); alIn = params(5); beIn = params(6); delIn = params(7); kapIn = params(8);
            OptionsIn = SnakeOptions.getSpecifiedOptions(DEBUG, iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn);
            
            iterOut = params(9); WlOut = params(10); WeOut = params(11); WtOut = params(12); alOut = params(13); beOut = params(14); delOut = params(15); kapOut = params(16);
            OptionsOut = SnakeOptions.getSpecifiedOptions(DEBUG, iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut);
        end
        
        function [OptionsIn, OptionsOut] = getLiverOptions(DEBUG)
            params = [1.6645   -0.6115    3.8542    0.0160    0.3467    0.0264   -0.1131    1.0516    3.8982    0.0516   -6.8314    0.3190   -0.0825    0.4198    0.6603   -2.2152]; % multiple objectives with 1-dice
%             params = [1.4276   -0.3007    2.4772    0.0108    0.2288    0.2311   -0.0613    1.7769    3.5745    0.0341    1.8300    0.1016    0.2066    0.2252    1.0686   -1.3575]; % multiple objectives with 1/dice and max
%             params = [1.4527   -0.4702    1.9121    0.0076    0.2302    0.2271   -0.0794    2.5498    3.9983    0.0258    2.7338    0.0626    0.1784    0.1640    0.8441   -1.2624];  % max penalty
%             params = [1.8582   -0.5917    1.7512    0.0101    0.2161    0.2068   -0.0473    2.0106    3.3873    0.0437    2.4288    0.0999    0.1396    0.1691    0.9379   -2.2118]; % max penalty no dice
%             init_it_in = 2; init_wline_in = -0.5; init_wedge_in = wedge_default; init_wterm_in = wterm_default; init_alpha_in = alpha_default;
%             init_beta_in = beta_default; init_delta_in = -0.05; init_kappa_in = kappa_default; 
%             init_it_out = 3; init_wline_out = wline_default; init_wedge_out = wedge_default; init_wterm_out = 0.1; init_alpha_out = alpha_default;
%             init_beta_out = beta_default; init_delta_out = 1; init_kappa_out = -2; 
%             params = [init_it_in, init_wline_in, init_wedge_in, init_wterm_in, init_alpha_in, init_beta_in, init_delta_in, init_kappa_in, ...
%                 init_it_out, init_wline_out, init_wedge_out, init_wterm_out, init_alpha_out, init_beta_out, init_delta_out, init_kappa_out]; % original parameters
            iterIn = params(1); WlIn = params(2); WeIn = params(3); WtIn = params(4); alIn = params(5); beIn = params(6); delIn = params(7); kapIn = params(8);
            OptionsIn = SnakeOptions.getSpecifiedOptions(DEBUG, iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn);
            
            iterOut = params(9); WlOut = params(10); WeOut = params(11); WtOut = params(12); alOut = params(13); beOut = params(14); delOut = params(15); kapOut = params(16);
            OptionsOut = SnakeOptions.getSpecifiedOptions(DEBUG, iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut);
        end
        
        function [OptionsIn, OptionsOut] = getKidneyOptions(DEBUG)
            OptionsIn = struct;
            OptionsIn.Iterations = 2;
            OptionsIn.Delta = -0.5;

            OptionsOut = struct;
            OptionsOut.Iterations = 2;
            OptionsOut.Delta = 1;
            OptionsOut.Kappa = -2;
            OptionsOut.Wterm = 0.1;
            
            if(DEBUG)
                OptionsIn.Verbose = true;
                OptionsOut.Verbose = true;
            end
        end
        
        function [OptionsIn, OptionsOut] = getBrainOptions(DEBUG)
            OptionsIn = struct;
            OptionsIn.Iterations = 5;
            OptionsIn.Delta = -0.5;

            OptionsOut = struct;
            OptionsOut.Iterations = 2;
            OptionsOut.Delta = 1;
            OptionsOut.Kappa = -2;
            OptionsOut.Wterm = 0.1;
            
            if(DEBUG)
                OptionsIn.Verbose = true;
                OptionsOut.Verbose = true;
            end
        end
        
        function [OptionsStruct] = getAllOptions(DEBUG)
            OptionsStruct = struct;
            
            % lung
            [OptionsIn, OptionsOut] = SnakeOptions.getLungOptions(DEBUG);
            OptionsStruct.LungOptions = {OptionsIn, OptionsOut};
            % liver
            [OptionsIn, OptionsOut] = SnakeOptions.getLiverOptions(DEBUG);
            OptionsStruct.LiverOptions = {OptionsIn, OptionsOut};
            % kidney
            [OptionsIn, OptionsOut] = SnakeOptions.getKidneyOptions(DEBUG);
            OptionsStruct.KidneyOptions = {OptionsIn, OptionsOut};
            % brain
            [OptionsIn, OptionsOut] = SnakeOptions.getBrainOptions(DEBUG);
            OptionsStruct.BrainOptions = {OptionsIn, OptionsOut};
        end
    end
    
end

