classdef Prior
    %PRIOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        %TODO: index of disparity for estimating if contour is good or not.
        % mean of contour grayscae divided by std
        function globalUncertainty = getGlobalUncertainty(seg,certaintyPrior, contourOnly)
            if ~exist('contourOnly','var')
                contourOnly = false;
            end
            
            if contourOnly
                segContour = logical(seg - imerode(seg,ones(3,3,3)));
                certaintyVals = certaintyPrior(segContour);
            else
                certaintyVals = certaintyPrior(seg > 0);
            end
            
            if sum(sum(sum(certaintyPrior>1)))>0 ||sum(sum(sum(certaintyPrior<0)))>0
                error 'illegal certainty map';
            end
            
            %high values = certainty. taht's why we need to flip the
            %results
            globalUncertainty = 1-mean(certaintyVals);
        end
        
        
        function priorsResCell = calculateAlgoPriors(im,algosCell,params)
            priorsResCell = cell(size(algosCell));
            for ii=1:length(algosCell)
                priorsResCell{ii} = Prior.calculatePriors(im,algosCell{ii},params);
            end
        end
        
        
        
        
        function priorRes = calculatePriors(im,seg,params)
            relevanceCell = cell(0,0);
            priorRes = AlgoPriorStrct;
            
            
            if isstruct(im)
                imNii = im;
                segNii = seg;
                im = imNii.img;
                seg = segNii.img;
            end
            priorRes.boundary = Utils.getBoundries(seg);
            if isfield(params,'intensity') && params.intensity
                if ~exist('params.filterSize','var')
                    priorRes.intensity = Prior.intensityPrior(im,seg);
                else
                    priorRes.intensity = Prior.intensityPrior(im,seg,params.filterSize);
                end
                relevanceCell = [relevanceCell priorRes.intensity];
            end
            if isfield(params,'intensityLocal') && params.intensityLocal
                if ~isfield(params,'intensityLocalT')
                    params.intensityLocalT = VariabilityExperiment.calculateThreshold(seg,im);
                end
                priorRes.intensityLocal = Prior.intensityPriorLocal(im,seg,params.filterSize,params.intensityLocalT);
                relevanceCell = [relevanceCell priorRes.intensityLocal];
            end
            if isfield(params,'curvature') && params.curvature
                priorRes.curvature = Prior.curvaturePrior(seg);
                relevanceCell = [relevanceCell priorRes.curvature];
            end
            if isfield(params,'shape') && params.shape
                if ~isfield(params,'shapePriorMask')
                    priorRes.shapePriorMask = Prior.registerAtlas(imNii,segNii,params.shapePriorAtlas);
                else
                    priorRes.shapePriorMask = params.shapePriorMask;
                end
                
                priorRes.shape = Prior.shapePriorFromRegisteredAtlas(seg,priorRes.shapePriorMask);
                %priorRes.shape = Prior.shapePriorGlobalHeuristics(seg);
                relevanceCell = [relevanceCell priorRes.shape];
            end
            if isfield(params,'texture') && params.texture
                priorRes.texture = Prior.texturePrior(im,seg);
                relevanceCell = [relevanceCell priorRes.texture];
            end
            if isfield(params,'min') && params.min
                priorRes.min = min(cat(4,relevanceCell{:}),[],4);
            end
            
            if isfield(params,'max') && params.max
                priorRes.max = max(cat(4,relevanceCell{:}),[],4);
            end
            if isfield(params,'mean') && params.mean
                priorRes.mean = mean(cat(4,relevanceCell{:}),[],4);
            end
        end
        
        function [res] = shapePriorFromRegisteredAtlas(seg,shapePriorMask)
            res = zeros(size(shapePriorMask));
            CLOSE_DIST_THRESH = 4;
            HIGH_DIST_THRESH = 8;
            distsMap = Prior.getDistsMap(shapePriorMask,CLOSE_DIST_THRESH,HIGH_DIST_THRESH);
            segContours = Utils.getBoundries(seg);
            for z=1:size(distsMap,3)
                currentRes = res(:,:,z);
                currentDist = distsMap(:,:,z);
                currentRes(segContours(:,:,z)) = currentDist(segContours(:,:,z));
                res(:,:,z) = currentRes;
            end
        end
        
        function [shapePriorMask] = registerAtlas(imNii,segNii,shapeAtlas)
            NII_FILE_NAME = 'temp\imNiiFileForEval.nii.gz';
            SEG_FILE_NAME = 'temp\segNiiFileForEval.nii.gz';
            
            if exist(NII_FILE_NAME,'file')==2
                delete(NII_FILE_NAME);
            end
            if exist(SEG_FILE_NAME,'file')==2
                delete(SEG_FILE_NAME);
            end
            
            if ischar(imNii)
                copyfile(imNii,NII_FILE_NAME);
            else
                imNii.untouch=0;
                save_nii_gzip(imNii,NII_FILE_NAME);
            end
            if ischar(segNii)
                copyfile(segNii,SEG_FILE_NAME);
                segNii = IO.loadFile(segNii);
            else
                segNii.untouch=0;
                save_nii_gzip(segNii,SEG_FILE_NAME);
            end
            
            [shapePriorMask, success] = shapeAtlas.registerToNii(NII_FILE_NAME, SEG_FILE_NAME);
            if ~success
                error('couldnt register atlas to nii');
            end
            
            delete(NII_FILE_NAME);
            delete(SEG_FILE_NAME);
        end
        
        function [mapCell, fullMapPriorsCell, algos] = evaluateAlgos(fileStrct, nAlgos, w)
            if ~exist('w','var')
                w = ones(4,1);
            end
            mapCell = cell(length(fileStrct.algoFiles),1);
            fullMapPriorsCell = cell(length(fileStrct.algoFiles),1);
            algos = cell(length(fileStrct.algoFiles),1);
            
            im = load_untouch_nii_gzip(fileStrct.file);
            im = im.img;
            if ~exist('nAlgos','var')
                nAlgos = length(fileStrct.algoFiles)
            end
            validInds = [];
            for i=1:min(length(fileStrct.algoFiles),nAlgos)
                [i length(fileStrct.algoFiles)]
                try
                    seg = load_untouch_nii_gzip(fileStrct.algoFiles{i});
                catch exception
                    continue;
                end
                seg = seg.img;
                [map, priorsCell] = Prior.getTotalPrior(im, seg, w);
                mapCell{i} = map;
                fullMapPriorsCell{i} = priorsCell;
                algos{i} = seg;
                validInds = [validInds i];
            end
            fullMapPriorsCell = fullMapPriorsCell(validInds);
            mapCell = mapCell(validInds);
            algos = algos(validInds);
            
        end
        
        function [map, priorsCell] = getTotalPrior(im, seg, w)
            intensityPrior = zeros(size(im));
            texturePrior = zeros(size(im));
            curvaturePrior = zeros(size(im));
            shapePrior = zeros(size(im));
            
            priorsCell = cell(1,1);
            if w(1) > 0
                intensityPrior = Prior.intensityPrior(im, seg);
                priorsCell{length(priorsCell)+1} = intensityPrior;
            end
            if w(2) > 0
                texturePrior = Prior.texturePrior(im, seg);
                priorsCell{length(priorsCell)+1} = texturePrior;
            end
            if w(3) > 0
                curvaturePrior = Prior.curvaturePrior(seg);
                priorsCell{length(priorsCell)+1} = curvaturePrior;
            end
            if w(4) > 0
                shapePrior = Prior.shapePrior(seg);
                priorsCell{length(priorsCell)+1} = shapePrior;
            end
            if ~exist('w','var')
                w = ones(3,1);
            end
            map = (w(1)*intensityPrior + w(2)*texturePrior + w(3)*curvaturePrior + w(4)*shapePrior) / sum(w);
            priorsCell{1} = map;
        end
        
        
        function res  = calcAngleInArrAux(xy,gap,ind)
            if ind<=gap || ind > size(xy,1)-gap
                res = -1;
            else
                res = VariabilityEstimator.calcAngle(xy(ind-gap,1),xy(ind-gap,2),xy(ind+gap,1),xy(ind+gap,2))+(pi/2);
            end
        end
        
        function outIm = borderRep(im,len)
            [m,n,z] = size(im);
            outIm = zeros(m+len*2,n+len*2,z);
            
            outIm(len+1:end-len,len+1:end-len,:) = im;
            
            outIm(1:len,len+1:end-len,:) = repmat(im(1,:,:),len,1);
            outIm(end-len+1:end,len+1:end-len,:) = repmat(im(end,:,:),len,1);
            
            outIm(len+1:end-len,1:len,:) = repmat(im(:,1,:),1,len);
            outIm(len+1:end-len,end-len+1:end,:) = repmat(im(:,end,:),1,len);
            
        end
        function [map] = intensityPriorLocal(im, seg,filterSize,T)
            if ~exist('filterSize','var')
                filterSize = 3;
            end
            gap = filterSize;
            im = Prior.borderRep(im,gap);
            im = mat2gray(im);
            seg = Prior.borderRep(seg,gap)>0;
            map = zeros(size(im));
            %mag = Utils.gradientMagDiffSizes(im,filterSize);
            boundries = Utils.getBoundries(seg);
            kernel = Utils.getSobelFilt(gap*2+1);
            kernel = kernel(gap+1,:);
            %testMask = zeros(size(boundries));
            for z=find(sum(sum(seg,1),2)>0)'
                currentIm = im(:,:,z);
                currentMap = zeros(size(currentIm));
                CC = bwconncomp(boundries(:,:,z),8);
                for tt=1:length(CC)
                    if length(CC.PixelIdxList{tt})<10
                        continue;
                    end
                    [y,x] = ind2sub(size(boundries(:,:,z)),CC.PixelIdxList{tt});
                    [yx] = bwtraceboundary(boundries(:,:,z),[y(1) x(1)],'W',8);
                    if isequal(yx(end,:),yx(1,:))
                        yx = yx(1:end-1,:);
                    end
                    
                    curvGap = 4;
                    yx = [yx(end-curvGap+1:end,:) ; yx ; yx(1:curvGap,:)];
                    y = yx(:,1);
                    x = yx(:,2);
                    
                    f = @(ind)VariabilityEstimator.calcAngle(yx(ind-curvGap,2),yx(ind-curvGap,1),yx(ind+curvGap,2),yx(ind+curvGap,1))+(pi/2);
                    angles = arrayfun(f,curvGap+1:size(yx,1)-curvGap);
                    for t=1:length(angles)

                        pntsArr1 = VariabilityEstimator.pointPathFromAngle(y(t),x(t),angles(t),gap);
                        pntsArr2 = VariabilityEstimator.pointPathFromAngle(y(t),x(t),angles(t)+pi,gap);
                        if size(pntsArr2,1)~=gap+1 || size(pntsArr1,1)~=gap+1 || ...
                                ~isequal(pntsArr1(gap+1,:),[y(t),x(t)]) ||  ~isequal(pntsArr2(gap+1,:),[y(t),x(t)])
                            error 'handle this';
                        end
                        pntsArr = [pntsArr1(1:end-1,:);pntsArr2(end:-1:1,:)];
                        pntsVals = currentIm(sub2ind(size(currentIm),pntsArr(:,1),pntsArr(:,2)));
                        currentMap(sub2ind(size(currentIm),y(t),x(t))) = abs(dot(pntsVals,kernel))/(sum(abs(kernel)/2));
                        %testMask = testMask | drawLines(size(testMask),pntsArr1(1,:),pntsArr1(end,:));
                        %testMask = testMask | drawLines(size(testMask),pntsArr2(1,:),pntsArr2(end,:));
                    end
                    %uncertaintyRegions.angle{tt} = VariabilityEstimator.calcAngle(x(1),y(1),x(end),y(end))+(pi/2);
                    
                    if ~exist('T') || isempty(T)
                        morphSize = 3;
                        dilatedBoun = Utils.getBoundries(imdilate(seg(:,:,z),strel('disk',morphSize)));
                        erodedSeg = imerode(seg,strel('disk',morphSize));
                        %erodedBoun = Utils.getBoundries(erodedSeg);
                        val1 = mean(currentIm(dilatedBoun));
                        val2 = mean(currentIm(erodedSeg));
                        T = abs(val1-val2)/2;
                    end
                    currentMap(currentMap>T) = 1;
                    map(:,:,z) = mat2gray(currentMap);
                end
            end
            


           
            %testMask = testMask*0.5;
            %testMask(boundries) = 1;
            map = map(gap+1:end-gap,gap+1:end-gap,:);
        end
        %intensity prior, where the thrshold is global
        function map = intensityPrior(im, seg,filterSize)
            
            if ~exist('filterSize','var')
                filterSize = 3;
            end
            
            map = zeros(size(im));
            mag = Utils.gradientMagDiffSizes(im,filterSize);
            
            prevMinEdgeVals = ones(size(im,3),1)*(-1);
            for i=1:size(im,3)
                currentMap = zeros(size(map(:,:,i)));
                
                if sum(sum(seg(:,:,i))) > 0
                    currentMag = mag(:,:,i);
                    currentIm = im(:,:,i);
                    
                    [Y, X] = ind2sub(size(currentMag),find(seg(:,:,i)));
                    
                    PADDING = 5;
                    minY = max(min(Y) - PADDING,1);
                    minX = max(min(X) - PADDING,1);
                    maxY = min(max(Y) + PADDING, size(seg,1));
                    maxX = min(max(X) + PADDING, size(seg,2));
                    clear X Y;
                    
                    currentMagOnRoi  = currentMag(minY:maxY,minX:maxX);
                    edges = edge(currentIm(minY:maxY,minX:maxX)) & currentMagOnRoi > 0.00001;
                    if sum(edges(:)) > 0
                        minEdgeVal = min(min(currentMagOnRoi(edges)));
                    else
                        minEdgeVal = max(currentMagOnRoi(:)) + 1;
                    end
                    
                    %a = sort(currentMagOnRoi(:));
                    
                    %minEdgeVal = a(int16(length(a)*0.95));
                    
                    prevMinEdgeVals(i) = minEdgeVal;
                    if i>2 &&prevMinEdgeVals(i-1)>0 &&prevMinEdgeVals(i-2)>0
                        minEdgeVal = min([minEdgeVal,prevMinEdgeVals(i-1),prevMinEdgeVals(i-2)]);
                    end
                    
                    
                    currentMagOnRoi(currentMagOnRoi > minEdgeVal) = minEdgeVal;
                    currentMagOnRoi = currentMagOnRoi ./ minEdgeVal;
                    currentMag(minY:maxY,minX:maxX) = currentMagOnRoi;
                    
                    contour = Utils.getBoundries(seg(:,:,i));
                    
                    %nonContour = ~contour & seg(:,:,i);
                    currentMap(contour) = currentMag(contour);
                    %currentMap = imdilate(currentMap,ones(3,3));
                    
                    %nonContorMag = zeros(size(currentMag));
                    %nonContorMag(nonContour) = 1-currentMag(nonContour);
                    %currentMap(nonContour) = nonContorMag(nonContour);
                end
                %currentMap = imdilate(currentMap,ones(3,3));
                map(:,:,i) = currentMap;
                
                
            end
            
            
            map(~seg) = 0;
        end
        
        
        
        function map = texturePrior2(im, seg)
            map = zeros(size(im));
            [Y, X, Z] = ind2sub(size(im),find(seg));
            filtR=generateRadialFilterLBP(8, 1);
            PADDING = 600;
            minY = max(min(Y) - PADDING,1);
            minX = max(min(X) - PADDING,1);
            minZ = max(min(Z) - PADDING,1);
            maxY = min(max(Y) + PADDING, size(seg,1));
            maxX = min(max(X) + PADDING, size(seg,2));
            maxZ = min(max(Z) + PADDING, size(seg,3));
            clear X Y Z;
            
            imInRoi = im(minY:maxY,minX:maxX,minZ:maxZ);
            segInRoi = seg(minY:maxY,minX:maxX,minZ:maxZ);
            mapInRoi = map(minY:maxY,minX:maxX,minZ:maxZ);
            for z=minZ:maxZ
                l = z - minZ + 1;
                currentIm = imInRoi(:,:,l);
                currentSeg = segInRoi(:,:,l);
                contour = Utils.getBoundries(currentSeg);
                
                currentMap= efficientLBP(currentIm, 'filtR', filtR);
                
                mapInRoi(:,:,l) = currentMap;
            end
            map(minY:maxY,minX:maxX,minZ:maxZ) = mapInRoi;
        end
        
        %regions with std above max will get 0. and below min will get 1.
        function map = texturePrior(im, seg, minStdThresh, maxStdThresh, contourMode)
            %if contourMode is off - the prior of contours will be 1
            if ~exist('contourMode', 'var')
                contourMode = true;
            end
            map = zeros(size(im));
            [Y, X, Z] = ind2sub(size(im),find(seg));
            
            PADDING = 5;
            minY = max(min(Y) - PADDING,1);
            minX = max(min(X) - PADDING,1);
            minZ = max(min(Z) - PADDING,1);
            maxY = min(max(Y) + PADDING, size(seg,1));
            maxX = min(max(X) + PADDING, size(seg,2));
            maxZ = min(max(Z) + PADDING, size(seg,3));
            clear X Y Z;
            
            imInRoi = im(minY:maxY,minX:maxX,minZ:maxZ);
            segInRoi = seg(minY:maxY,minX:maxX,minZ:maxZ);
            
            %std filters on 3x3 env, of each layer
            stdVals = stdfilt(imInRoi);
            stdVals(~segInRoi) = 0;
            
            if ~exist('minStdThresh','var')
                
                KIDNEY_MIN_VAR_THRESHOLD = 35;
                KIDNEY_MAX_VAR_THRESHOLD = 100;
                minStdThresh = KIDNEY_MIN_VAR_THRESHOLD;
                maxStdThresh = KIDNEY_MAX_VAR_THRESHOLD;
            end
            stdVals = min(stdVals,maxStdThresh);
            stdVals = max(stdVals,minStdThresh);
            stdVals = (stdVals - minStdThresh)/(maxStdThresh-minStdThresh);
            
            map(minY:maxY,minX:maxX,minZ:maxZ) = stdVals;
            map = 1-map;
            
            %sets contour prior to be 1.
            
            if ~contourMode
                se = zeros(3,3,3);
                se(:,:,2) = strel('disk',3);
                contour = logical(seg - imerode(seg,se));
                map(contour) = 1;
            else
                
                for z=minZ:maxZ
                    currentMap = map(:,:,z);
                    currentSeg = seg(:,:,z);
                    
                    contour = Utils.getBoundries(currentSeg);
                    distFromContor = bwdist(contour,'chessboard');
                    d = 5;
                    relevanceMask = distFromContor >=d & distFromContor <=(d+1) & currentSeg;
                    contourRes = ones(size(distFromContor));
                    contourRes(relevanceMask) = currentMap(relevanceMask);
                    se = ones(d*2+1,d*2+1);
                    contourRes = imerode(contourRes,se);
                    currentMap(distFromContor<3) = contourRes(distFromContor<3);
                    
                    
                    currentMap(~currentSeg) = 0;
                    currentMap(~contour) = 1;
                    map(:,:,z) = currentMap;
                end
                
            end
            %noise reduction
            %map = imdilate(map,se);
            map(~seg) = 0;
        end
        
        function contDist = getContourDist(layer, threshLow, threshHigh)
            contour = Utils.getBoundries(layer);
            
            contDist = bwdist(contour);
            contDist(contDist > threshHigh) = threshHigh;
            contDist(contDist < threshLow) = threshLow;
            contDist = contDist - threshLow;
            contDist = 1-contDist/(threshHigh-threshLow);
            
        end
        
        
        function map = shapePriorGlobalHeuristics(seg,minThresh,maxThresh)
            
            if ~exist('minThresh','var')
                minThresh = 5;
            end
            if ~exist('maxThresh','var')
                maxThresh = 15;
            end
            
            map = zeros(size(seg));
            SMALL_REGION_THRESH = 50;
            relevantLayers = find(sum(sum(seg,1),2) > SMALL_REGION_THRESH)';
            for z = relevantLayers
                %first step - finds the largest connected component
                currentSeg = seg(:,:,z);
                CC = bwconncomp(currentSeg,4);
                numPixels = cellfun(@numel,CC.PixelIdxList);
                [~,biggestCCInd] = max(numPixels);
                currentMap = zeros(size(currentSeg));
                %second stage - calculate the distance from the contour of
                %the biggest CC covex hull
                biggestCCMask = zeros(size(currentSeg));
                biggestCCMask(CC.PixelIdxList{biggestCCInd}) = 1;
                currentBoundry = Utils.getBoundries(biggestCCMask);
                currentMap(currentBoundry) = 1;
                K = bwconvhull(currentBoundry);
                contourPrior = Prior.getContourDist(K, minThresh, maxThresh);
                currentMap(currentBoundry) = contourPrior(currentBoundry);
                
                
                map(:,:,z) = currentMap;
                
            end
            
            map2 = Prior.shapePriorLayersConsistensy(seg);
            map = min(map,map2);
            
        end
        
        
        function distsMap = getDistsMap(seg,CLOSE_DIST_THRESH,HIGH_DIST_THRESH,relevantLayers)
            
            if ~exist('relevantLayers','var')
                relevantLayers = squeeze(sum(sum(seg,1),2) > 0);
            end
            layersToCalculate = unique(find(relevantLayers));
            
            distsMap = zeros(size(seg));
            contoursMap = Utils.getBoundries(seg);
            for z=layersToCalculate
                distsMap(:,:,z) = bwdist(contoursMap(:,:,z),'chessboard');
            end
            distsMap(distsMap > HIGH_DIST_THRESH) = HIGH_DIST_THRESH;
            distsMap(distsMap < CLOSE_DIST_THRESH) = CLOSE_DIST_THRESH;
            distsMap = distsMap - CLOSE_DIST_THRESH;
            distsMap = 1-distsMap/(HIGH_DIST_THRESH-CLOSE_DIST_THRESH);
        end
        
        function map = shapePriorLayersConsistensy(seg, contourOnly)
            if ~exist('contourOnly')
                contourOnly = true;
            end
            if contourOnly == false
                error 'mode not implemented';
            end
            
            CLOSE_DIST_THRESH = 4;%4^0.5;
            HIGH_DIST_THRESH = 8;%36^0.5;
            SMALL_REGION_THRESH = 50;
            map = zeros(size(seg));
            
            relevantLayers = find(sum(sum(seg,1),2) > SMALL_REGION_THRESH)';
            %filters out irrelevant indices
            deleteList = [];
            c = 1;
            for z=relevantLayers
                if sum(ismember(relevantLayers,z-1)) == 0 || sum(ismember(relevantLayers,z+1))==0
                    deleteList = [deleteList c];
                end
                c = c + 1;
            end
            
            relevantLayers(deleteList) = [];
            
            distsMap = Prior.getDistsMap(seg,CLOSE_DIST_THRESH,HIGH_DIST_THRESH,relevantLayers);
            contoursMap = Utils.getBoundries(seg);
            for z = relevantLayers
                %z
                prevLayerDist = distsMap(:,:,z-1);
                nextLayerDist = distsMap(:,:,z+1);
                combDist = (prevLayerDist + nextLayerDist)/2;
                currentContour =  contoursMap(:,:,z);
                currnetMapResult = zeros(size(currentContour));
                currnetMapResult(currentContour) = combDist(currentContour);
                map(:,:,z) = currnetMapResult;
            end
            
            if(max(map(:))) > 1 || min(map(:))<0
                error('invalid map');
            end
            
        end
        
        function [xOut, yOut, xy] = orderXY(x,y)
            xy = [x,y];
            dists = zeros(length(x)-1,1);
            for i=2:length(x)
                nextXY = xy(i-1,:);
                nextRep = repmat(nextXY,length(x) - i + 1,1);
                distsMap = abs(nextRep - xy(i:end,:));
                distsMap = sum(distsMap,2);
                minInd = find(distsMap==min(distsMap));
                minInd = minInd(1) + (i-1);
                xy([i,minInd(1)],:) = xy([minInd(1),i],:);
                dists(i-1) = min(distsMap);
            end
            
            
            xOut = xy(1,:);
            yOut = xy(2,:);
        end
        
        
        function map = curvaturePrior(seg, halfWinSize)
            if ~exist('halfWinSize','var')
                halfWinSize = 7;
            end
            
            %boundries = bwboundaries(seg);
            map = zeros(size(seg));
            THRESHOLD_MIN = 0.3;
            THRESHOLD_MAX = 2;
            
            for i=1:size(seg,3)
                currentMap = zeros(size(map(:,:,i)));
                currentSeg = seg(:,:,i);
                
                if sum(currentSeg(:)) > 0
                    %inner segmentation part gets highest prior
                    currentMap(currentSeg > 0) = 0;
                    
                    %calculate curvature on boundries
                    [curvatureArr,xy] = Utils.calcCurvatore(currentSeg,halfWinSize,THRESHOLD_MIN,THRESHOLD_MAX);
                    
                    %noise cleaning
                    curvatureArr = imdilate(curvatureArr,ones(5,1));
                    
                    %updates map
                    currentMap(sub2ind(size(currentSeg),xy(:,1),xy(:,2))) = 1-curvatureArr;
                    
                end
                map(:,:,i) = currentMap;
                
            end
            
            %boundries = seg - imerode(seg,ones(3,3,3));
            
            %[Y, X, Z] = ind2sub(size(boundries),find(boundries));
            %[normals,curvature] = findPointNormals(points,[],[0,0,10],true);
            %[K,H,Pmax,Pmin] = surfature(X,Y,Z);
            
            
            %map = 1;
        end
        function [mag, dx, dy] = imgradient (im)
            %function [xderiv,yderiv] = imgradient (im, sigma)
            %gk = gaussKernel(sigma);
            %dg = dgausskernel(sigma);
            %padsize = floor(length(gk)/2);
            %padim = padarray(im, [padsize padsize], 'replicate');
            
            %xderiv = conv2(gk, dg, double(padim), 'valid');
            %yderiv = conv2(dg, gk, double(padim), 'valid');
            
            [dx, dy] = gradient(double(im));
            mag = (dx.^2 + dy.^2).^0.5;
        end
        
    end
    
end