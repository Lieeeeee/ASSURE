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
        function [OptionsIn, OptionsOut] = getLungOptions(DEBUG)
            OptionsIn = struct;
            OptionsIn.Iterations = 5;
%             OptionsIn.Delta = 1;
%             OptionsIn.Wline = -0.04;
%             OptionsIn.Kappa = -2;

            OptionsOut = struct;
            OptionsOut.Iterations = 5;
            OptionsOut.Delta = 1;
            OptionsOut.Wline = -0.04;
            OptionsOut.Kappa = -2;
            
            if(DEBUG)
                OptionsIn.Verbose = true;
                OptionsOut.Verbose = true;
            end
        end
    end
    
end
