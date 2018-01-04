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
            params = [2.5451    0.0334    2.1619    0.0088    0.1600    0.1711   -0.5740    2.3761    3.5504   -0.0586   -2.0904   -0.0071    0.1349    0.2608    0.4343    3.0997];
            
            iterIn = params(1); WlIn = params(2); WeIn = params(3); WtIn = params(4); alIn = params(5); beIn = params(6); delIn = params(7); kapIn = params(8);
            OptionsIn = SnakeOptions.getSpecifiedOptions(DEBUG, iterIn, WlIn, WeIn, WtIn, alIn, beIn, delIn, kapIn);
            
            iterOut = params(9); WlOut = params(10); WeOut = params(11); WtOut = params(12); alOut = params(13); beOut = params(14); delOut = params(15); kapOut = params(16);
            OptionsOut = SnakeOptions.getSpecifiedOptions(DEBUG, iterOut, WlOut, WeOut, WtOut, alOut, beOut, delOut, kapOut);
        end
        
        function [OptionsIn, OptionsOut] = getLiverOptions(DEBUG)
            OptionsIn = struct;
            OptionsIn.Iterations = 2;
            OptionsIn.Delta = -0.05;
            OptionsIn.Wline = -0.5;

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

