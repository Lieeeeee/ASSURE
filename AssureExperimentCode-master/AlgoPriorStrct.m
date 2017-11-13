classdef AlgoPriorStrct
    %ALGOPRIOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        intensity
        intensityLocal
        curvature
        shape
        texture
        min
        max
        mean
        boundary
        mainPrior
        shapePriorMask
    end
    methods(Static)
        
        
        function params = getParams()
            params.intensity = true;
            params.curvature = true;
            params.shape = false;
            params.texture = true;
            params.max = false;
            params.min = true;
            params.mean = false;
        end
        
        function cell = setCellMainPrior(cell,mainPrior)
            for ii=1:length(cell)
                cell{ii}.mainPrior = mainPrior;
            end
        end
    end
    methods
        function strct = AlgoPriorStrct
            strct.mainPrior = 'min';
        end
        
        function [globalErr, globalErrPerLayer] = getGlobalErr(strct)
            prior = strct.getMainPrior();
            globalErr = mean(prior(strct.boundary));
            globalErrPerLayer = zeros(size(strct.boundary,3),1);
            for z=1:size(strct.boundary,3)
               currentLayer = prior(:,:,z);
               if sum(sum(strct.boundary(:,:,z)))==0
                  globalErrPerLayer(z) = nan;
                  continue;
               end
               globalErrPerLayer(z) =  mean(currentLayer(strct.boundary(:,:,z)));
            end
            
        end
        function prior = getMainPrior(strct)
            if ~any(strcmp(properties(strct), strct.mainPrior))
               error(['field: ' strct.mainPrior ' not exist']);
            else 
                prior = strct.(strct.mainPrior);
            end    
        end
        
        function prior = filteredMainPrior(strct)
            if ~any(strcmp(properties(strct), strct.mainPrior))
               error(['field: ' strct.mainPrior ' not exist']);
            else 
                prior = strct.(strct.mainPrior);
            end    
            MED_SIZE = 5;
            for z=find(sum(sum(strct.boundary,1),2)>0,1,'first'):find(sum(sum(strct.boundary,1),2)>0,1,'last')
                mask = strct.boundary(:,:,z);
                currPrior = prior(:,:,z);
                currPrior(~mask) = nan;
                currPrior = nanmedfilt2(currPrior,[MED_SIZE,MED_SIZE]);
                currPrior(isnan(currPrior)) = 0;
                prior(:,:,z)=currPrior;
            end

        end
    end
    
end

