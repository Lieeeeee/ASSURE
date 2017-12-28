classdef VariabilityEstimator
    %VARIABILITYESTIMATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f
        eps
        im
        seg
        uncertaintyT
        winLength
    end
    
    
    methods (Static)
        function [x1,y1] = pointFromAngle(x0,y0,angle,lineLength)
            x1 = int16(x0 + lineLength * cos(angle));
            y1 = int16(y0 + lineLength * sin(angle));
        end
        
        function pntsArr = pointPathFromAngle(y0,x0,angle,lineLength)
            done = false;
            modifiedLength = lineLength;
            c = 0;
            while ~done
                [x1,y1] = VariabilityEstimator.pointFromAngle(x0,y0,angle,modifiedLength);
                [pntsArr] = generateLinesPoints([y1,x1], [y0,x0]);
                
                actualLength = size(pntsArr,1)-1;
                done = actualLength == lineLength;
                if actualLength < lineLength
                    modifiedLength = modifiedLength + 1;
                elseif actualLength > lineLength
                    pntsArr = pntsArr(end-lineLength:end,:);
                    done = true;
                end
                c = c +1;
                if c >5
                    error 'infinite loop'
                end
            end
        end
        
        function [angle] = calcAngle(x0,y0,x1,y1)
            angle = atan((y1-y0)./(x1-x0));
        end
        
        function [newSeg,roi] = generateNewSeg(seg,indOut,indIn)
            roi = false(size(seg));
            roi(indOut) = false;
            roi(indIn) = true;
            newSeg = Utils.getBoundries(seg,false);
            newSeg(indOut) = 0;
            newSeg(indIn) = 1;
        end
        function [meanChange, mainPriorOut] = calculatePriorsWrapper(im,seg,params,mainPrior,roi)
            priorRes = Prior.calculatePriors(im,seg,params);
            priorRes.mainPrior = mainPrior;
            mainPriorOut = priorRes.getMainPrior();
            
            if exist('roi','var') && ~isempty(roi)
                meanChange = mean(mainPriorOut(roi));
            else
                meanChange = mean(mainPriorOut(seg));
            end
        end
        
        
        
        function [varMask, varData,sensitivityMask] = evaluate3DVarMaskTrivial(im,seg,priorRes,params,eps)
            varData = cell(size(im,3),1);
            varMask = false(size(seg));
            sensitivityMask = zeros([size(seg,1),size(seg,2),size(seg,3),3]);
            for z = find(sum(sum(seg,1),2) > 0)'
                varMask(:,:,z) = imdilate(seg(:,:,z),strel('disk',2)) &~ imerode(seg(:,:,z),strel('disk',2)); %TODO
            end
        end
        
        function [varMask, varData,sensitivityMask] = evaluate3DVarMask(im,seg,priorRes,params,eps)
            varData = cell(size(im,3),1);
            varMask = false(size(seg));
            sensitivityMask = zeros([size(seg,1),size(seg,2),size(seg,3),3]);
            for z = find(sum(sum(seg,1),2) > 0)'
               % z
                currentIm = im(:,:,z);
                currentSeg = seg(:,:,z);
                currentParams = params;
                if isfield(currentParams,'shapePriorMask')
                    currentParams.shapePriorMask = currentParams.shapePriorMask(:,:,z);
                end
                func = @(segmentationIn, roi) VariabilityEstimator.calculatePriorsWrapper(currentIm,segmentationIn,currentParams,'min',roi);
                if isfield(params,'winLength')
                    estimator = VariabilityEstimator(currentIm,currentSeg,func,eps,[],params.winLength);
                else
                    estimator = VariabilityEstimator(currentIm,currentSeg,func,eps,[]);
                end
                
                [uncertaintyRegions] = estimator.extractUncertaintyRegions(priorRes(:,:,z));
                [variabilityData] = estimator.measureRegionsSensitivity(uncertaintyRegions);
                
                if ~isfield(params,'minUncertaintyRegionSize')
                    params.minUncertaintyRegionSize = 0;
                end
                variabilityData = estimator.evaluateVariability(variabilityData,params.minUncertaintyRegionSize);
                currentVariabilityMask =  estimator.getVariabilityMask(variabilityData);
                varMask(:,:,z) = currentVariabilityMask;
                varData{z} = variabilityData;
                currentSensitivityMask = getSensitivityDbgRes(estimator, variabilityData);
                sensitivityMask(:,:,z,:) = currentSensitivityMask;
            end
        end
        
        function [varMask, varData, sensitivityMask] = evaluate3DVarMaskActiveContours(im, seg)
            % evaluates a variability mask, given a segmentation and its priors. 
            varData = cell(size(im,3),1);
            varMask = false(size(seg));
            sensitivityMask = zeros([size(seg,1),size(seg,2),size(seg,3),3]);
            
            alg = 'edge'; % 'Chan-Vese'; 
            n_iterations = 4; % 5
            contraction_param_out = -0.5; % -0.5
            smooth_param_out = 0.1; % 0
            contraction_param_in = 1;
            smooth_param_in = 1;
            for z = find(sum(sum(seg,1),2) > 0)'
                currentIm = im(:,:,z);
                currentSeg = seg(:,:,z);
                
                % first, apply active contour outwards, starting from the given segmentataion
                bw_1 = activecontour(currentIm, currentSeg, n_iterations, alg, 'ContractionBias', contraction_param_out, 'SmoothFactor', smooth_param_out);
                
                % next, apply active contour inwards, starting from the given segmentataion
                bw_2 = activecontour(currentIm, currentSeg, n_iterations, alg, 'ContractionBias', contraction_param_in, 'SmoothFactor', smooth_param_in);
                
                varMask(:,:,z) = (bw_1 | bw_2) - (bw_1 & bw_2); % bw_1 - bw_2
            end
        end
        
        function [varMask, varData, sensitivityMask] = evaluate3DVarMaskSnakes(im, seg, OptionsIn, OptionsOut)
            % evaluates a variability mask, given a segmentation and its priors. 
            varData = cell(size(im,3),1);
            sensitivityMask = zeros([size(seg,1),size(seg,2),size(seg,3),3]);
            varMask = false(size(seg));
            
            for z = find(sum(sum(seg,1),2) > 0)'
                currentSeg = seg(:,:,z);
                currentIm = im(:,:,z);
                
                % prepare the points
                [currentSegPoints] = Utils.getPointsOnContour(currentSeg);
                
                % update options.nPoints
%                 OptionsIn.nPoints = floor(size(currentSegPoints,1)/3);
                OptionsIn.nPoints = size(currentSegPoints,1);
                OptionsOut.nPoints = size(currentSegPoints,1);
                
                % first, apply active contour outwards, starting from the given segmentataion
                [~, bw_1] = Snake2D(currentIm,fliplr(currentSegPoints),OptionsOut);
                
                % next, apply active contour inwards, starting from the given segmentataion
                [~, bw_2] = Snake2D(currentIm,fliplr(currentSegPoints),OptionsIn);
                
                varMask(:,:,z) = (bw_1 | bw_2 | currentSeg) - (bw_1 & bw_2); % bw_1 - bw_2
            end
        end
        
        function out = calcMeanShape(seg1,seg2)
            out = false(size(seg1));
            
            for z = find(sum(sum(seg1&seg2,1),2) > 0)'
                s1 = seg1(:,:,z);
                s2 = seg2(:,:,z);
                inter = s1 & s2;
                CC = bwconncomp(s1);
                
                for pp=1:length(CC.PixelIdxList)
                    o = out(:,:,z);
                    [y,x] = ind2sub(size(inter),CC.PixelIdxList{pp}(1));
                    s1CC = Utils.getCurrentPixelCC(s1,[y,x]);
                    [Y1,X1] = ind2sub(size(s1CC),find(s1CC));
                    [yx1] = bwtraceboundary(s1,[Y1(1),X1(1)],'W',8);
                    
                    s1BndryMask = false(size(s1CC));
                    s1BndryMask(sub2ind(size(s1BndryMask),yx1(:,1),yx1(:,2))) = 1;
                    
                    s2Inter = s1CC & s2;
                    if sum(s2Inter(:))==0
                        continue;
                    end
                    [Y2,X2] = ind2sub(size(s2Inter),find(s2Inter));
                    s2CC = Utils.getCurrentPixelCC(s2,[Y2(1),X2(1)]);
                    [Y2,X2] = ind2sub(size(s2CC),find(s2CC));
                    [yx2] = bwtraceboundary(s2,[Y2(1),X2(1)],'W',8);
                    s2BndryMask = false(size(s2CC));
                    s2BndryMask(sub2ind(size(s2BndryMask),yx2(:,1),yx2(:,2))) = 1;
                    
                    % x=54 y=26
                    fNorm = @(p1,P) sum((repmat(p1,size(P,1),1)-P).^2,2);
                    fMinDist = @(p1,P)find(fNorm(p1,P)==min(fNorm(p1,P)),1,'first');
                    yxNew = [];
                    for tt=1:size(yx1,1)
                        p = yx1(tt,:);
                        cp = yx2(fMinDist(p,yx2),:);
                        avP = int16((p+cp)/2);
                        yxNew = [yxNew;avP];
                        
                    end
                    o(sub2ind(size(o),yxNew(:,1),yxNew(:,2)))=true;
                    
                    lineAdditionMask = zeros(size(o));
                    for tt=1:size(yxNew,1)-1
                        p1 = yxNew(tt,:);
                        p2 = yxNew(tt+1,:);
                        if max(abs(p1-p2))>=2
                            % lineAdditionMask = lineAdditionMask | drawLines(size(lineAdditionMask),p1(end:-1:1),p2(end:-1:1));
                            lineAdditionMask = lineAdditionMask | drawLines(size(lineAdditionMask),p1,p2);
                        end
                    end
                    o = o | lineAdditionMask;
                    o = imfill(o,'holes');
                    out(:,:,z) =  out(:,:,z) | o;
                end
                
                %closes holes in contour
                
                
            end
            
        end
        
        function out = calcMeanShapeMultSegs(segs)
            [union,intersection] = Utils.calcUnionIntersection(segs);
            out = VariabilityEstimator.calcMeanShape(union,intersection);
            
        end
    end
    
    methods
        
        
        %function grad = calcCurveGradient(yx)
        %    t=1;
        %    maxLength = 10;
        %    y = yx(t:5:end,1);
        %    x = yx(t:5:end,2);
        %    dx  = gradient(x);
        %    dy  = gradient(y);
        %    angles = atan(dy./dx)+(pi/2);
        %    grad = angles;
        %
        %    x2 = x - maxLength * cos(angles);
        %    y2 = y - maxLength * sin(angles);
        %
        %end
        function estimator = VariabilityEstimator(im, seg, estimationFunction, sensitivityVal,uncertaintyT,winLength)
            estimator.f = estimationFunction;
            estimator.eps = sensitivityVal;
            estimator.im = im;
            estimator.seg = seg;
            if ~exist('uncertaintyT','var') || isempty(uncertaintyT)
                uncertaintyT = 0.7;
            end
            
            if ~exist('winLength','var')
                winLength = 3;
            end
            estimator.uncertaintyT = uncertaintyT; %0.5;
            estimator.winLength = winLength;
        end
        
        function [res,sensitivityPnts] = measureRegionSensitivity(estimator,regionInds,angle)
            [y,x] = ind2sub(size(estimator.seg),regionInds);
            %TODO make sure that angle + winLength is inside window scope
            [x0,y0] = VariabilityEstimator.pointFromAngle(x,y,angle,estimator.winLength);
            [m, n] = size(estimator.seg);
            
            legitInd = find((x0<=n & y0<=m & x0>=1 & y0>=1),1,'last');
            x0 = x0(1:legitInd);
            y0 = y0(1:legitInd);
            
            [~,sensitivityPnts] = drawLines(size(estimator.seg),[y,x],[y0,x0]);
            res = zeros(size(sensitivityPnts,1),1);
            for rr=1:size(sensitivityPnts,1)
                xVals = sensitivityPnts(rr,2,:);
                yVals = sensitivityPnts(rr,1,:);
                [newSeg,roi] = VariabilityEstimator.generateNewSeg(estimator.seg,regionInds,sub2ind(size(estimator.seg),yVals,xVals));
                res(rr)=estimator.f(newSeg,roi);
            end
        end
        
        function variabilityData = measureRegionsSensitivity(estimator,uncertaintyRegions)
            nRegions = length(uncertaintyRegions.PixelIdxList);
            resIn = cell(nRegions,1);
            sensitivityPntsIn = cell(nRegions,1);
            resOut = cell(nRegions,1);
            sensitivityPntsOut = cell(nRegions,1);
            for tt=1:nRegions
                [resIn{tt},sensitivityPntsIn{tt}] = measureRegionSensitivity(estimator,uncertaintyRegions.PixelIdxList{tt},uncertaintyRegions.angle{tt});
                [resOut{tt},sensitivityPntsOut{tt}] = measureRegionSensitivity(estimator,uncertaintyRegions.PixelIdxList{tt},uncertaintyRegions.angle{tt}+pi);
            end
            variabilityData = VariabilityData;
            variabilityData.resIn = resIn;
            variabilityData.resOut = resOut;
            variabilityData.sensitivityPntsIn = sensitivityPntsIn;
            variabilityData.sensitivityPntsOut = sensitivityPntsOut;
            variabilityData.uncertaintyRegions = uncertaintyRegions;
        end
        
        function [uncertaintyRegions,lowPriorMask] = extractUncertaintyRegions(estimator,prior)
            UNCERTAINTY_T = estimator.uncertaintyT;
            CC_AREA_T = 3;
            
            boundries = Utils.getBoundries(estimator.seg);
            distMap = bwdist(boundries);
            distMap(~estimator.seg) = -inf;
            [Cy,Cx] = ind2sub(size(distMap),find(max(distMap(:))));
            
            %filter min
            prior(~boundries) = inf;
            prior = imerode(prior,strel('disk',2));
            prior(~boundries) = 0;
            %finds CC, seperate when there's high curvature
            [curv,xy] = Utils.calcCurvatore(estimator.seg, 1);
            xy(curv<1,:)=[];
            prior(sub2ind(size(prior),xy(:,1),xy(:,2))) = inf;
            lowPriorMask = prior < UNCERTAINTY_T & boundries;
            
            uncertaintyRegions = bwconncomp(lowPriorMask,8);
            rp = regionprops(uncertaintyRegions, 'Area', 'PixelIdxList');
            uncertaintyRegions.PixelIdxList([rp.Area]<CC_AREA_T)= [];
            
            for tt=1:length(uncertaintyRegions.PixelIdxList)
                
                [y,x] = ind2sub(size(prior),uncertaintyRegions.PixelIdxList{tt});
                [yx] = bwtraceboundary(lowPriorMask,[y(1) x(1)],'W',8);
                if isequal(yx(end,:),yx(1,:))
                    yx = yx(1:end-1,:);
                end
                y = yx(:,1);
                x = yx(:,2);
                
                uncertaintyRegions.angle{tt} = VariabilityEstimator.calcAngle(x(1),y(1),x(end),y(end))+(pi/2);
                
                %the angle should point to the outern part of the
                %segmentation
                uncertaintyRegions.cInd{tt} = int16(length(x))/2;
                angleToC = VariabilityEstimator.calcAngle(Cx,Cy,x(uncertaintyRegions.cInd{tt}),y(uncertaintyRegions.cInd{tt}));
                if abs(angleToC-uncertaintyRegions.angle{tt}) > pi/2
                    uncertaintyRegions.angle{tt} = uncertaintyRegions.angle{tt}-pi;
                end
                
            end
        end
        
        function ind = evaluateRegionVariability(estimator,res)
            res = res > estimator.uncertaintyT;
            firstCertain = find(res,1,'first');
            if isempty(firstCertain)
                ind = length(res);
            else
                firstNonCertein = find(~res(firstCertain+1:end),1,'first');
                if isempty(firstNonCertein)
                    ind = length(res);
                else
                    ind = firstNonCertein + firstCertain - 1; %-1 because we exclude the uncertainty part
                end
            end
        end
        
        function variabilityData = evaluateVariability(estimator,variabilityData,minRegionSize)
            if ~exist('minRegionSize','var')
               minRegionSize = 0; 
            end
            variabilityData.sensitivityEndIn = zeros(size(variabilityData.resIn));
            variabilityData.sensitivityEndOut = zeros(size(variabilityData.resOut));
            for tt=1:length(variabilityData.resIn)
                if length(variabilityData.uncertaintyRegions.PixelIdxList{tt})<minRegionSize
                    continue;
                end
                ind = estimator.evaluateRegionVariability(variabilityData.resIn{tt});
                variabilityData.sensitivityEndIn(tt) = ind;
            end
            for tt=1:length(variabilityData.resOut)
                if length(variabilityData.uncertaintyRegions.PixelIdxList{tt})<minRegionSize
                    continue;
                end
                ind = estimator.evaluateRegionVariability(variabilityData.resOut{tt});
                variabilityData.sensitivityEndOut(tt) = ind;
            end
            
            combRes = [variabilityData.sensitivityEndIn,variabilityData.sensitivityEndOut];
            combRes = min(combRes')';
            variabilityData.sensitivityEndOut = combRes;
            variabilityData.sensitivityEndIn = combRes;
            %CLOSE_DIST = 3; %TODO optimize
            %for tt=1:length(variabilityData.resOut)
            %    if variabilityData.sensitivityEndIn(tt) <= CLOSE_DIST && variabilityData.sensitivityEndOut(tt) > CLOSE_DIST
            %        variabilityData.sensitivityEndOut(tt) = variabilityData.sensitivityEndIn(tt);
            %    elseif variabilityData.sensitivityEndOut(tt) <= CLOSE_DIST && variabilityData.sensitivityEndIn(tt) > CLOSE_DIST
            %       variabilityData.sensitivityEndIn(tt) = variabilityData.sensitivityEndOut(tt);
            %    end
            %end
            
        end
        function variabilityMask = getVariabilityMask(estimator,variabilityData)
            variabilityMask = false(size(estimator.im));
            pntsAll = [variabilityData.sensitivityPntsIn;variabilityData.sensitivityPntsOut];
            stopPntsAll = [variabilityData.sensitivityEndIn;variabilityData.sensitivityEndOut];
            stopPntsAll = max(2,stopPntsAll);
            for tt=1:length(pntsAll)
                pnts = pntsAll{tt};
                %p1 = pnts(1,:,1);
                %p2 = pnts(1,:,size(pnts,3));
                %p3 = pnts(stopPntsAll(tt),:,1);
                %p4 = pnts(stopPntsAll(tt),:,size(pnts,3));
                %addMask = drawLines(size(variabilityMask),p1,p2);
                %addMask = addMask | drawLines(size(addMask),p1,p3);
                %addMask = addMask | drawLines(size(addMask),p2,p4);
                %addMask = addMask | drawLines(size(addMask),p3,p4);
                %variabilityMask = variabilityMask | bwfill(addMask,'holes');
                
                tempMask = zeros(size(variabilityMask));
                for zz=1:stopPntsAll(tt)
                    g = min(zz,estimator.winLength);
                    bLeft = g;
                    bRight = size(pnts,3)-g+1;
                    
                    tempMask(sub2ind(size(variabilityMask),squeeze(pnts(zz,1,bLeft:bRight)),squeeze(pnts(zz,2,bLeft:bRight)))) = 1;
                    %variabilityMask(sub2ind(size(variabilityMask),squeeze(pnts(zz,1,:)),squeeze(pnts(zz,2,:)))) = 1;
                end
                tempMask = imclose(tempMask,strel('disk',3));
                variabilityMask = variabilityMask | tempMask;
            end
            
            diffMask = imdilate(estimator.seg,strel('disk',2)) &~ imerode(estimator.seg,strel('disk',2)); %TODO
            variabilityMask = variabilityMask | diffMask;
            %variabilityMask = imclose(variabilityMask,strel('disk',3));
            %  variabilityMask = imopen(variabilityMask,strel('disk',1));
        end
        
        
        
        function res = getSensitivityDbgRes(estimator, variabilityData)
            
            
            resMask = zeros(size(estimator.im));
            changeMask = zeros(size(estimator.im));
            for tt=1:length(variabilityData.uncertaintyRegions.PixelIdxList)
                
                cInd = variabilityData.uncertaintyRegions.cInd{tt};
                centralLineYX = variabilityData.sensitivityPntsIn{tt}(:,:,cInd);
                %lineMask = zeros(size(resMask));
                %lineMask(sub2ind(size(lineMask),centralLineYX(:,1),centralLineYX(:,2))) = 1;
                
                resMask(sub2ind(size(resMask),centralLineYX(:,1),centralLineYX(:,2))) = variabilityData.resIn{tt};
                changeMask(sub2ind(size(resMask),centralLineYX(:,1),centralLineYX(:,2))) = 1;
                
                centralLineYX = variabilityData.sensitivityPntsOut{tt}(:,:,cInd);
                
                resMask(sub2ind(size(resMask),centralLineYX(:,1),centralLineYX(:,2))) = variabilityData.resOut{tt};
                changeMask(sub2ind(size(resMask),centralLineYX(:,1),centralLineYX(:,2))) = 1;
            end
            res = Utils.displayCertaintyUncertainty2_3D(estimator.im,resMask,Utils.getBoundries(changeMask,false));
            res = Utils.displaySegmentation(res,estimator.seg,[0,1,1]);
            
            %figure,IO.imshow3D(res);
            
        end
        
        function displayUncertaintyRegions(estimator,priorIm,variabilityData)
            %imshow(prior); hold on;
            if isempty(variabilityData.sensitivityEndOut)
                sensitivityEndOut = cellfun(@length,variabilityData.resOut);
                sensitivityEndIn = cellfun(@length,variabilityData.resIn);
            else
                sensitivityEndOut = variabilityData.sensitivityEndOut;
                sensitivityEndIn = variabilityData.sensitivityEndIn;
            end
            
            
            figure,imshow(priorIm); hold on;
            for tt=1:length(variabilityData.uncertaintyRegions.PixelIdxList)
                list = variabilityData.uncertaintyRegions.PixelIdxList{tt};
                angle = variabilityData.uncertaintyRegions.angle{tt};
                [Cy,Cx] = ind2sub(size(priorIm),list(int16(length(list))/2));
                [x0,y0] = VariabilityEstimator.pointFromAngle(Cx,Cy,angle,sensitivityEndIn(tt));
                [x1,y1] = VariabilityEstimator.pointFromAngle(Cx,Cy,angle+pi,sensitivityEndOut(tt));
                
                quiver(Cx,Cy,Cx-x0,Cy-y0,'Color','blue','lineWidth',1.5,'MaxHeadSize',50);
                quiver(Cx,Cy,Cx-x1,Cy-y1,'Color','red','lineWidth',1.5,'MaxHeadSize',50);
                %lowPriorMask = lowPriorMask | drawLines(size(lowPriorMask),[y0,x0],[y1,x1]);
            end
            
        end
        
        %CURRENTLY NOT IN USED
        function [contourLines, simpilfyContour] = getContourPatches(estimator, seg)
            nLines = 50;
            [Y,X] = ind2sub(size(seg),find(seg));
            [yx] = bwtraceboundary(seg,[Y(1) X(1)],'W',8);
            anchors = simplifyPoly(yx, nLines);
            
            % determines the points of each line, and it's boundries
            [simpilfyContour, linesPnts] = drawLines(size(seg),anchors,[anchors(end,:);anchors(1:end-1,:)]);
            contourLines.pnts = linesPnts;
            contourLines.bndryPnts = cell(nLines,1);
            contourLines.bndryPnts{1} = anchors;
            contourLines.bndryPnts{2} = [anchors(end,:);anchors(1:end-1,:)];
            
            %calculates curvature for each line
            yx = [anchors(end,:);anchors;anchors(1,:)];
            dx  = gradient(yx(:,2));
            dy  = gradient(yx(:,1));
            contourLines.dirs = atan(dy./dx) +(pi/2);
            contourLines.dirs = contourLines.dirs(2:end-1);
            
            %determines centreal point for each line
            contourLines.Cx = contourLines.pnts(length(contourLines.pnts),2);
            contourLines.Cy = contourLines.pnts(length(contourLines.pnts),1);
            
            contourLines.x0 = contourLines.Cx + maxLength * cos(angles);
            contourLines.y0 = contourLines.Cy + maxLength * sin(angles);
            contourLines.x1 = contourLines.Cx - maxLength * cos(angles);
            contourLines.y1 = contourLines.Cy - maxLength * sin(angles);
            %simplfyContour = zeros(size(seg));
            
            
        end
        
        function estimateLineVariability(estimator,seg,lib)
            
            
        end
        
        function res = estimate(estimator, im, seg)
            contour = Utils.getBoundries(seg);
            res = contour;
            maxLength = 10;
            for z = find(sum(sum(seg,1),2)>0,1,'first'):find(sum(sum(seg,1),2)>0,1,'last')
                currSeg = seg(:,:,z);
                currentIm = im(:,:,z);
                currSeg = Utils.keepBiggestCC(currSeg);
                
                [contourLines, simpilfiedContour] = getContourPatches(estimator, currSeg);
                
                imshow(simpilfiedContour);
                
            end
            
        end
        
        
        
        %  function res = estimate(estimator, im, seg)
        %      %TODO - pad!
        %      contour = Utils.getBoundries(seg);
        %      mag = Utils.calcGradientMag(im,seg);
        %      res = contour;
        %      halfWinSize = 3;
        %      maxLength = 10;
        %      for z = find(sum(sum(seg,1),2)>0,1,'first'):find(sum(sum(seg,1),2)>0,1,'last')
        %          currSeg = seg(:,:,z);
        %          currentIm = im(:,:,z);
        %          currSeg = Utils.keepBiggestCC(currSeg);
        %
        %          [Y,X] = ind2sub(size(currSeg),find(currSeg));
        %          [xy] = bwtraceboundary(currSeg,[Y(1) X(1)],'W',8);
        %          for t=1:halfWinSize
        %              y = xy(t:halfWinSize:end,1);
        %              x = xy(t:halfWinSize:end,2);
        %              dx  = gradient(x);
        %              dy  = gradient(y);
        %              angles = atan(dy./dx) +(pi/2);
        %
        %              x0 = x + maxLength * cos(angles);
        %              y0 = y + maxLength * sin(angles);
        %              x1 = x - maxLength * cos(angles);
        %              y1 = y - maxLength * sin(angles);
        %          end
        
        
        %     end
        
        % end
        

    end
    
end



