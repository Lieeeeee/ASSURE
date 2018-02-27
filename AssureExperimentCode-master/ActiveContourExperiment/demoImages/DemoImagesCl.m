classdef DemoImagesCl
    %DEMOIMAGESCL class for creating and handling synthetic images for
    %testing
    %   Detailed explanation goes here
    
    properties (Constant)
        circleRadius = 15;
        circleLocation = [50 50];
    end
    
    methods (Static)
        function [I] = getDarkCircleImg()
            I = zeros(100,100);
            % add circle in center of image
            I = rgb2gray(insertShape(I, 'FilledCircle', [DemoImagesCl.circleLocation(1) DemoImagesCl.circleLocation(2) DemoImagesCl.circleRadius], 'Color', 'white', 'Opacity', 1));
        end
        
        function [I] = getLightCircleImg()
            I = DemoImagesCl.getDarkCircleImg();
            I = 1-I;
        end
        
        function [blurred_I] = addBlurToImg(I, sigma)
            if ~exist('sigma','var')
                sigma = 3;
            end
            blurred_I = imgaussfilt(I, sigma);
        end
        
        function [seg] = getCircleSeg(radius, location)
            seg = zeros(100,100);
            seg = rgb2gray(insertShape(seg, 'FilledCircle', [location(1) location(2) radius], 'Color', 'white', 'Opacity', 1));
            seg = seg > 0.5;
        end
        
        function [I, seg] = getCircleAndSeg(isDark)
            if ~exist('isDark','var')
                isDark = true;
            end
            if isDark
                I = DemoImagesCl.getDarkCircleImg();
            else
                I = DemoImagesCl.getLightCircleImg();
            end
            I = DemoImagesCl.addBlurToImg(I);
            seg = DemoImagesCl.getCircleSeg(DemoImagesCl.circleRadius, DemoImagesCl.circleLocation);
        end
    end
    
    methods

    end
    
end

