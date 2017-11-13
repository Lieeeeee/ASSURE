classdef Unitest
    %UNITEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        function out = jaccard(A,B)
           andMask = A & B;
           out = sum(andMask(:)) / (sum(A(:)) + sum(B(:)) - sum(andMask(:)));
        end
        
        function out = dice(A,B)
           andMask = A & B;
           out = (2*sum(andMask(:))) / (sum(A(:)) + sum(B(:)));
        end

        function cellOut = sortAlgoCellByJaccardRes(algoCell,gt)
            errStrct = Unitest.calculateJaccardRes(algoCell,gt);
            [~,I] = sort(errStrct.jaccard3D);
            cellOut = algoCell(I);
        end
        
        function errStrct = calculateGlobalErrPrior(algoStrct,errStrct)
           if ~iscell(algoStrct)
               error 'input should be cell'
           end
           if ~exist('errStrct','var')
               errStrct = [];
           end
           
           errStrct.globalErr = zeros(size(algoStrct));
           errStrct.globalErrPerLayer = cell(size(algoStrct));
           
           for ii=1:length(algoStrct)
               [globalErr, globalErrPerLayer] = getGlocalErr(algoStrct{ii});
                errStrct.globalErr(ii) = globalErr;
                errStrct.globalErrPerLayer{ii} = globalErrPerLayer;
           end
           
        end
        
        function errStrct = contourMatchingErr(algo,gt,radius)
            if ~exist('radius','var')
                radius = 15;
            end
            
            s1 = Utils.getBoundries(algo);
            s2 = Utils.getBoundries(gt);
            
            [Y1,X1] = find(s1);
            s1Points = [Y1,X1];
            [Y2,X2] = find(s2);
            s2Points = [Y2,X2];
            
            for t = 1:length(s1Points)
               p = s1Points(t,:);
                
            end
            
        end
        
        function errStrct = calculateJaccardRes(algos,gt,errStrct)
           if ~iscell(algos)
               error 'input should be cell'
           end
           if ~exist('errStrct','var')
               errStrct = [];
           end
           
           errStrct.jaccard2D = cell(size(algos));
           errStrct.jaccard3D = zeros(size(algos));
           
           for ii=1:length(algos)
              errStrct.jaccard3D(ii) = Unitest.jaccard(algos{ii},gt);
              errStrct.jaccard2D{ii} = zeros(size(algos{ii},3),1);
              for tt = 1:size(algos{ii},3)
                  errStrct.jaccard2D{ii}(tt) = Unitest.jaccard(algos{ii}(:,:,tt),gt(:,:,tt));
              end
           end
           
        end
        function CorrelationStrct = getCorelationResults(file, nAlgos)
            gtSegFile = file.seg;
            %reads ground truth
            gtSeg = load_untouch_nii_gzip(gtSegFile); gtSeg = gtSeg.img;
            origNii = load_untouch_nii_gzip(file.file);
            
            
            %my algo
            [~,fullMapPriorsCell,algos] = Prior.evaluateAlgos(file,nAlgos);
            nPriors = length(fullMapPriorsCell{1});
            nAlgos = length(algos);
            %
            %gt eval
            [uncertaintyDistGTCell, intensityGTCell, segSampleCell, jackardArr] = Unitest.calcUncertaintyGTCell(algos, gtSeg, origNii.img);
            jackardArr
            
            %for each algo, calculates the correlation between each prior
            %and GT method
            localCorrelationTable = zeros(nPriors,2,nAlgos);
            
            %for each prior calculoates the correlation between actual
            %error and estimated error
            %                      surfaceDist, intensityDist,
            %
            % allPriors
            % intensityPrior
            % texturePrior
            %curvaturePrior
            globalCorrelationTable = zeros(nPriors,1);
            
            
            for p=1:nPriors
                globalErrors = zeros(size(jackardArr));
                for j=1:nAlgos
                    [p  nPriors j nAlgos]
                    %algo segmentation
                    algoSeg = algos{j};
                    %uncertainty estimation
                    certaintyMap = fullMapPriorsCell{j}{p};
                    %global error
                    globalErrors(j) = Prior.getGlobalUncertainty(algoSeg,certaintyMap,true);
                    

                    
                    if max(certaintyMap(:)) > 1
                       error 'illegal map prior'; 
                    end
                   
                    %extracts local uncertainty values
                    uncertaintyVals = 1 - certaintyMap(segSampleCell{j});
                    
                    %ground truth estimation values
                    intensityGTVals = intensityGTCell{j} / max(intensityGTCell{j});
                    distGTVals = uncertaintyDistGTCell{j} / max(uncertaintyDistGTCell{j});
                    
                    localCorrelationTable(p,1,j) = corr(uncertaintyVals,distGTVals);
                    localCorrelationTable(p,2,j) = corr(uncertaintyVals,intensityGTVals);
                    
                end
                globalCorrelationTable(p) = corr(globalErrors,1-jackardArr);
            end
            localCorrelationTable(isnan(localCorrelationTable)) = 0;
            CorrelationStrct.globalCorrelationTable = globalCorrelationTable;
            CorrelationStrct.localCorrelationTable = localCorrelationTable;
            errorCorr = corr(globalErrors,1-jackardArr);
            CorrelationStrct.errorCorr = errorCorr;
            
            
            
        end
        
        
        
        function [surfDistCell, intensityDistCell, segSampleCell, jackardArr] = calcUncertaintyGTCell(segCell, gtSeg, img)
            surfDistCell = cell(length(segCell),1);
            intensityDistCell = cell(length(segCell),1);
            segSampleCell = cell(length(segCell),1);
            jackardArr = zeros(length(segCell),1);
            
            for i=1:length(segCell)
                if isempty(segCell{i})
                    break;
                end
                segSampleMask = Unitest.sampleSegmentationPoints(segCell{i},gtSeg);
                resSD = Unitest.calcSurfaceDistFromGT(segSampleMask, gtSeg);
                if sum(isnan(resSD)) > 0
                    a = 5;
                end
                surfDistCell{i} = resSD;
                segSampleCell{i} = segSampleMask;
                intensityDistCell{i} = Unitest.calcIntensityDiffFromGT(segSampleMask, gtSeg,img );
                andMask = gtSeg & segCell{i};
                jackardArr(i) = sum(andMask(:)) / (sum(gtSeg(:)) + sum(segCell{i}(:)) - sum(andMask(:)));
            end
            
            %segSampleCell = segSampleCell(1:i-1);
            %intensityDistCell = intensityDistCell(1:i-1);
            %surfDistCell = surfDistCell(1:i-1);
            %jackardArr = jackardArr(1:i-1);
            
        end
        
        
        
        function [segSampleMask] = sampleSegmentationPoints(seg, gtSeg, jumpFactor)
            if ~exist('jumpFactor','var')
                jumpFactor = 5;
            end
            seg = seg>0;
            segSampleMask = zeros(size(seg));
            relevantLayers = find(sum(sum(seg,1),2) > 0 & sum(sum(gtSeg,1),2) > 0)';
            layersToIgnore = [];
            c = 0;
            for z=relevantLayers
                c = c + 1;
                curGt = gtSeg(:,:,z);
                curSeg = seg(:,:,z);
                andMask = curGt & curSeg;
                jackardRes= sum(andMask(:)) / (sum(curGt(:)) + sum(curSeg(:)) - sum(andMask(:)));
                if jackardRes < 0.5
                    %imshow(mat2gray(curGt));figure,imshow(mat2gray(curSeg));
                    %figure,imshow(mat2gray(curSeg));
                    layersToIgnore = [layersToIgnore c];
                    continue;
                end
                CCSeg = bwconncomp(curSeg);
                CCGt = bwconncomp(curGt);
                if length(CCSeg.PixelIdxList) ~= length(CCGt.PixelIdxList)
                    firstCCSegRatio = length(CCSeg.PixelIdxList{1}) / sum(curSeg(:));
                    if (firstCCSegRatio < 0.99 &&  length(CCGt.PixelIdxList)==1) ||  length(CCGt.PixelIdxList)>1 
                        layersToIgnore = [layersToIgnore c];
                    end
                    continue;
                end
            end
            layersToIgnore
            relevantLayers(layersToIgnore) = [];
            
            %TODO: ignoring gt layers
            for z=relevantLayers
                
                currentSeg = seg(:,:,z);
                pointsStrct = bwboundaries(currentSeg,8,'noholes');
                %TODO - pointsStrct{1} should be pointsStrct{t}
                for t=1:length(pointsStrct)
                    
                    points = pointsStrct{1};
                    points = points(1:jumpFactor:end,:);
                    ind = sub2ind(size(currentSeg),points(:,1),points(:,2));
                    mask = zeros(size(currentSeg));
                    mask(ind) = 1;
                    segSampleMask(:,:,z) = mask;
                end
                
            end
            segSampleMask = logical(segSampleMask);
        end
        
        
        function outPoints = findClosestPoints(P1,P2)
            % inputs: P1 and P2 are lists of 2D points of sizes mx2 and nx2.
            % output: outPoints - matrix of size mx2, s.t outPoints(:,
            IDX = knnsearch(P2,P1,'Distance','euclidean','K',1);
            outPoints = P2(IDX,:);
            
        end
        
        function displayLayersAndErrs(combIm,errStrct,layer)
           [m,n,~,~] = size(combIm);
           nPics = length(errStrct.jaccard2D);
           y = ones(nPics,1)*(m - 20);
           x = ([((n/nPics)/4):(n/nPics):n]);
          
           figure,imshow(squeeze(combIm(:,:,layer,:)));
           for tt=1:nPics
              str = ['J: ' num2str(errStrct.jaccard2D{tt}(layer),2)];
              str = [str ',   G:' num2str(errStrct.globalErrPerLayer{tt}(layer),2)];
              
              text(x(tt),y(tt),str,'Color',[1 1 1],'FontSize',18);
           end

        end
        
        function [res, P1, P2] = calcIntensityDiffFromGT(samplesSegMask, gtSeg,volume)
            
            gtContour = logical(gtSeg - imerode(gtSeg,ones(3,3,3)));
            
            P1 = ones(sum(sum(sum(samplesSegMask>0))),3)*(-1);
            P2 = ones(size(P1))*(-1);
            res = ones(size(P1,1),1)*(-1);
            
            relevantLayers = find(sum(sum(samplesSegMask,1),2) > 0)';
            for z=relevantLayers
                segLayer = samplesSegMask(:,:,z);
                gtContourLayer = gtContour(:,:,z);
                if sum(sum(segLayer>0))<15 ||   sum(sum(gtContourLayer>0))<15
                    continue;
                end
                
                currentVolLayer = volume(:,:,z);
                
                [segY, segX] = ind2sub(size(segLayer),find(segLayer));
                segPoints = [segY,segX];
                [gtY, gtX] = ind2sub(size(gtContourLayer),find(gtContourLayer));
                
                gtPoints = [gtY,gtX];
                clear segX segY gtX gtY;
                closestPointsToSeg = Unitest.findClosestPoints(segPoints,gtPoints);
                
                nPoints = size(segPoints,1);
                firstNonInit = find(P1(:,1)==-1,1,'first');
                P1(firstNonInit:firstNonInit+nPoints-1,:) = [segPoints,ones(nPoints,1)*z];
                P2(firstNonInit:firstNonInit+nPoints-1,:) = [closestPointsToSeg,ones(nPoints,1)*z];
                
                P1Vals = currentVolLayer(sub2ind(size(currentVolLayer),segPoints(:,1),segPoints(:,2)));
                P2Vals = currentVolLayer(sub2ind(size(currentVolLayer),closestPointsToSeg(:,1),closestPointsToSeg(:,2)));
                res(firstNonInit:firstNonInit+nPoints-1) = (P1Vals-P2Vals).^2;
            end
            
        end
        
        function res = calcSurfaceDistFromGT(samplesSegMask, gtSeg)
            
            
            res = ones(sum(sum(sum(samplesSegMask>0))),1)*(-1);
            gtContour = logical(gtSeg - imerode(gtSeg,ones(3,3,3)));
            
            
            nLayers = sum(sum(sum(gtSeg,2),1)>0);
            maxDist = (sum(gtSeg(:))/nLayers)^0.5 / 2;
            relevantLayers = find(sum(sum(samplesSegMask,1),2) > 0)';
            
            for z=relevantLayers
                currentGT = gtContour(:,:,z);
                segLayer = samplesSegMask(:,:,z);
                
                distFromGT = bwdist(currentGT);
                distFromGT(distFromGT > maxDist) = maxDist;
                %distFromGT = max(distFromGT(:)) - distFromGT;
                distFromGT = distFromGT / max(distFromGT(:));
                
                resOnThisLayer = distFromGT(segLayer);
                nPoints = length(resOnThisLayer);
                firstNonInit = find(res==-1,1,'first');
                res(firstNonInit:firstNonInit + nPoints - 1) = resOnThisLayer;
            end
            
            
            
            
        end
        
    end
    
end

