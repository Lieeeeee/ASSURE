classdef VariabilityExperiment
    %VARIABILITYEXPERIMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function mapObj = getRadiologistDistributionMap(masks,masksFileNames,maxNAnnotators)
            sumMask = Utils.calcSegCellSum(masks);
            strCell = cell(size(masks{1}));
            for k=1:prod(size(strCell))
                strCell{k} = {};
            end
            
            for k=1:length(masks)
                [~,tempName,~] = fileparts(masksFileNames{k});
                name = SUREExperimentCls.getRadiologistName(tempName);
                curRelevanceMask = masks{k} & (sumMask <= maxNAnnotators);
                for kk=find(curRelevanceMask)'
                    strCell{kk} = {strCell{kk}{:},name};
                end
            end
            
            mapObj = containers.Map();
            for ind=find(sumMask<=maxNAnnotators & sumMask>0)'
                if isempty(strCell{ind})
                    error 'cell shouldnt be empty';
                else
                    str = strjoin(sort(strCell{ind}),',');
                    if mapObj.isKey(str)
                        mapObj(str) = mapObj(str) + 1;
                    else
                        mapObj(str) = 1;
                    end
                end
                
            end
        end
        
        function mapObj = getVarDistributionMap(masks,nAnnotatorsArr)
            sumMask = Utils.calcSegCellSum(masks);
            
            
            mapObj = containers.Map();
            for k=1:length(nAnnotatorsArr)
                nAnnotators = nAnnotatorsArr(k);
                str = num2str(nAnnotators);
                mapObj(str) = sum(sumMask(:)==nAnnotators);
            end
            
            
            %totalVoxels = sum(relevanceMask(:)>0);
            %for k=keys(mapObj)
            %    key = k{1};
            %    if isequal(key(1:2),', ')
            %        temp = (mapObj(key));
            %        mapObj.remove(key);
            %        key = key(3:end);
            %        mapObj(key) = temp;
            %    end
            %    mapObj(key) = mapObj(key) / totalVoxels;
            %end
            %vals = cell2mat(values(mapObj));
            %if abs(sum(vals)-1)>0.0000000001
            %    error 'sum should be 1'
            %end
        end
        
        function displayMapAsPie(mapObj,percentThrshold)
            %sortedVals = sort(cell2mat(values(mapObj)));
            %threshold = sortedVals(find(cumsum(sortedVals)>0.02,1,'first'));
            tempSum = 0;
            for k=keys(mapObj)
                if mapObj(k{1})<percentThrshold
                    tempSum = tempSum + mapObj(k{1});
                    mapObj.remove(k{1});
                end
            end
            if tempSum > 0
                mapObj('others') = tempSum;
            end
            
            vals = cell2mat(values(mapObj));
            %if abs(sum(vals)-1)>0.0000000001
            %    error 'sum should be 1'
            %end
            labels = keys(mapObj);
            labels = strcat(labels,':');
            
            h = pie(vals);
            hText = findobj(h,'Type','text'); % text object handles
            percentValues = get(hText,'String'); % percent values
            combinedtxt = strcat(labels',percentValues);
            combinedtxt(vals<0.005) = {''};
            for z=1:length(hText)
                hText(z).String = combinedtxt(z);
            end
            
        end
        
        function [m,n] = getSubPlotDim(nPlots)
            if nPlots >6
                m = 3; n = 3;
            elseif nPlots >4
                m = 2; n = 3;
            else
                m = 1; n = nPlots;
            end
        end
        
        function avgMap = getMapAvg(imgStrcts, mapObjCell)
            indsToIgnore = [];
            for z=2:length(imgStrcts)
                if length(imgStrcts{z}.masks) ~= length(imgStrcts{1}.masks)
                    str = ['example ' num2str(z) 'is out. diff num of annotations'];
                    warning(str);
                    indsToIgnore = [indsToIgnore,z];
                end
            end
            imgStrcts(indsToIgnore) = [];
            mapObjCell(indsToIgnore) = [];
            %get all keys
            mapLabels = {};
            for k=1:length(mapObjCell)
                curKeys = keys(mapObjCell{k});
                mapLabels = {mapLabels{:}, curKeys{:}};
            end
            [~,X,~] = unique(mapLabels,'stable');
            mapLabels = mapLabels(X);
            
            avgMap = containers.Map();
            for label=mapLabels
                curSum = 0;
                N = length(imgStrcts);
                label = label{1};
                for k=1:length(imgStrcts)
                    if mapObjCell{k}.isKey(label)
                        curSum = curSum + mapObjCell{k}(label);
                    end
                end
                avgMap(label) = curSum/N;
            end
        end
        
        function map = concatMaps(mapA,mapB)
            map = containers.Map(keys(mapA),values(mapA));
            for key=keys(mapB)
                key = key{1};
                if map.isKey(key)
                    map(key) = map(key) +  mapB(key);
                else
                    map(key) = mapB(key);
                end
            end
        end
        
        function mapOut = normalizeMap(mapIn)
            mapOut = containers.Map(keys(mapIn),values(mapIn));
            mapSum = sum(cell2mat(values(mapOut)));
            for key=keys(mapOut)
                key = key{1};
                mapOut(key) = mapOut(key)/mapSum;
            end
            
            if abs(sum(cell2mat(values(mapOut)))-1)>0.0000000001
                error 'sum should be 1'
            end
        end
        
        
        function avgMap = plotVarDistibution(imgStrcts,plotRadiologists, titleStr, plotConsensus )
            if ~exist('plotRadiologists','var')
                plotRadiologists = true;
            end
            if ~exist('plotConsensus','var')
                plotConsensus = false;
            end
            if ~exist('titleStr','var')
                titleStr = '';
            end
            if ~plotConsensus
                mainTitle = ['distribution of number of annotators dissagree with majority within variability - ' titleStr];
            else
                mainTitle = ['distribution of number of annotators dissagree with majority within possible - ' titleStr];
            end
            
            [m,n] = VariabilityExperiment.getSubPlotDim(length(imgStrcts)+1);
            
            mapsCell = {};
            for z=1:length(imgStrcts)
                subplot(m,n,z);
                sumMask = Utils.calcSegCellSum(imgStrcts{z}.masks);
                nAnnotators = length(imgStrcts{z}.masks);
                if plotConsensus
                    mapObj = VariabilityExperiment.getVarDistributionMap(imgStrcts{z}.masks,1:nAnnotators);
                    expectedSum = sum(sumMask(:)>0 & sumMask(:) <= nAnnotators );
                else
                    if plotRadiologists
                        initialN = 2;
                    else
                        initialN = 1;
                    end
                    mapObj = VariabilityExperiment.getVarDistributionMap(imgStrcts{z}.masks,initialN:ceil((nAnnotators-1)/2));
                    
                    contrastMasks = imgStrcts{z}.masks;
                    for zz=1:length(contrastMasks)
                        contrastMasks{zz} = ~contrastMasks{zz} & (sumMask > 0);
                    end
                    contrastMapObj = VariabilityExperiment.getVarDistributionMap(contrastMasks,initialN:floor((nAnnotators-1)/2));
                    mapObj = VariabilityExperiment.concatMaps(mapObj,contrastMapObj);
                    if plotRadiologists
                        mapRadObj = VariabilityExperiment.getRadiologistDistributionMap(imgStrcts{z}.masks,imgStrcts{z}.masksFilesNames,1);
                        contrastMapRadObj = VariabilityExperiment.getRadiologistDistributionMap(contrastMasks,imgStrcts{z}.masksFilesNames,1);
                        mapObj = VariabilityExperiment.concatMaps(mapObj,mapRadObj);
                        mapObj = VariabilityExperiment.concatMaps(mapObj,contrastMapRadObj);
                    end
                    
                    
                    expectedSum = sum(sumMask(:)>0 & sumMask(:) < nAnnotators );
                end
                mapObjSum = sum(cell2mat(values(mapObj)));
                if abs(mapObjSum - expectedSum) >0.0000000001
                    error 'expected sum is wrong'
                end
                
                
                mapObj =  VariabilityExperiment.normalizeMap(mapObj);
                
                mapsCell{z} = containers.Map(keys(mapObj),values(mapObj));
                VariabilityExperiment.displayMapAsBar(mapObj);
                xlabel('number of annotators');
                ylabel('percentage');
                %VariabilityExperiment.displayMapAsPie(mapObj,percentThrshold);
                [~,file,~] = fileparts(imgStrcts{z}.imgFileName);
                fileName = file(1:min(length(file),10));
                title(['file name: ' strrep(fileName,'_','-')]);
            end
            
            
            subplot(m,n,z+1);
            avgMap = VariabilityExperiment.getMapAvg(imgStrcts, mapsCell);
            %if abs(sum(cell2mat(values(avgMap)))-1)>0.0000000001
            %    error 'sum should be 1'
            %end
            VariabilityExperiment.displayMapAsBar(avgMap);
            title(['average - ' titleStr])
            
            annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, ...
                'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
            % bar( cell2mat( values(avgMap) ) )
            % set(gca,'XTick',[1:length(keys(avgMap))])
            % set(gca,'xticklabel', keys(avgMap))
            
        end
        
        function seniorityMap = getSeniorityMap(labels)
            isNumber = @(x)all(ismember(x, '0123456789+-.eEdD'));
             seniorityMap = containers.Map();
             for z=1:length(labels)
                if ~isNumber(labels{z})
                    seniorityMap(labels{z}) = SUREExperimentCls.getRadiologistSeniority(labels{z});
                else
                    error 'radiologist not found'
                end
             end

        end
        
        function displayMapAsBar(mapObj, groupIndsMap, colorByGroup, sortByGroup, groupNameStrct, linesMap)


            labels = keys(mapObj);
            vals = cell2mat(values(mapObj));
            vals = round(vals,2);
            
            if exist('linesMap','var')
               lineMapVals =  values(linesMap);
            end
            if ~exist('colorByGroup','var')
                colorByGroup = false;
            end
            if ~exist('sortByGroup','var')
                sortByGroup = false;
            end
            
            if exist('groupIndsMap','var')
               groupInds = cell2mat(values( groupIndsMap));
            end
            
            if exist('groupInds','var') && length(unique(groupInds))<=5
                colors = [1,0,0; 0,1,0; 0,0,1; 0,1,1; 1,1,0];
            else
                colors = Utils.getRGBColors();
            end
          
            
            if sortByGroup
                sortParam = groupInds;
                %colors = colors(sortParam+2,:);
                [~,I] = sort(sortParam);
            else
                I = 1:length(vals);
                %[~,I] = sort(vals);
                
                %[~,INames] = sort(labels);
                %colors = colors(INames,:);
            end
            
            if colorByGroup
                colors = colors(groupInds,:);
            elseif exist('groupColorVec','var')
                colors = colors(groupColorVec,:);
            else
                [~,INames] = sort(labels);
                colors = colors(INames,:);
            end
            
            
            hold on;
            
            
            %colors = {'r','g','b','y','c'};
            for z=1:length(I)
                ind = I(z);
                h=bar(z,vals(ind));
                if exist('linesMap','var')
                    plot([z-0.2,z+0.2],[lineMapVals{ind}(1) lineMapVals{ind}(1)],'Color',[0.5,0.5,0.5]);
                    plot([z-0.2,z+0.2],[lineMapVals{ind}(2) lineMapVals{ind}(2)],'Color',[0.5,0.5,0.5]);
                    plot([z,z],lineMapVals{ind}, 'Color',[0.7,0.7,0.7]);
                end
                set(h,'FaceColor',colors(ind,:));
                text(z,vals(ind),num2str(round(vals(ind),2)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9)
            end
            set(gca, 'XTick', 1:length(I), 'XTickLabel', labels(I));
            
            if exist('groupInds','var')
                
                sortedGroups = sort(groupInds);
                curGroupStart = 1;
                curGroupEnd = 1;
                curGroup = sortedGroups(curGroupStart);
                
                
                legendStrings  = {};
                h = zeros(length(unique(groupInds)), 1); c = 0;
                for z=(curGroupStart+1):length(sortedGroups)+1
                    if z<=length(sortedGroups) && sortedGroups(z) == curGroup
                        curGroupEnd = curGroupEnd + 1;
                    else
                        c = c +1 ;
                        currentMin = min(vals(I(curGroupStart:curGroupEnd)));
                        currentMax = max(vals(I(curGroupStart:curGroupEnd)));
                        currentMean = mean(vals(I(curGroupStart:curGroupEnd)));
                        curStrMain = ['mean:' num2str(round(currentMean,2)) '.'];
                        curCol = colors(I(curGroupStart),:);
                        %text(mean([curSeniorityStart,curSeniorityEnd]),currentMax+0.08,curStr,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'color',curCol)
                        curStrGap = ['-' num2str(round(currentMean-currentMin,2)) ',+' num2str(round(currentMax-currentMean,2))];
                        
                        if exist('groupNameStrct','var')
                            prefixStr = groupNameStrct{c};
                        else
                            prefixStr = '';
                        end
                        legendStrings = {legendStrings{:},[prefixStr ': ' curStrMain  curStrGap]};
                        
                        %text(mean([curSeniorityStart,curSeniorityEnd]),currentMax+0.06,curStr,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10,'color',curCol)
                        h(c) = plot(0,0,'or', 'visible', 'off','Color',curCol);
                        if z<=length(groupInds)
                            curGroupStart = z;
                            curGroupEnd = z;
                            curGroup = sortedGroups(z);
                            
                        end
                    end
                end
                legend(h,legendStrings);
                %legend(legendStrings);
            end
            
            if min(vals) >= 0 && max(vals)<=1
                if max(vals) > 0.5 && min(vals<0.8)
                    ylim([0 1]);
                elseif max(vals) > 0.5 && min(vals>=0.8)
                    ylim([0.8 1]);
                elseif max(vals)>0.3
                    ylim([0 0.5]);
                else
                    ylim([0 0.3]);
                end
            end
            hold off;
        end
        
        
        function avgAgr = getAvgAgreement(mask,others)
            resVec = zeros(length(others),1);
            for z=1:length(others)
                resVec(z) = Unitest.dice(mask,others{z});
            end
            avgAgr = mean(resVec);
            
        end
        
        function [mapAvg groupIndsMap] = rankCasesVar(imgCell, groupNames)
            if ~exist('measure','var')
                measure = 'diceFromOthers';
            end
            if ~iscell(imgCell{1})
                error 'expected for cell of cells';
            end

            m=3;n=3;
            mapObj = containers.Map();
            groupIndsMap =  containers.Map();
            mapCell = {};
            
            for tt=1:length(imgCell)
                for rr=1:length(imgCell{tt})
                    rr
                    
                    meanVal = 0;
                    for k=1:length(imgCell{tt}{rr}.masks)
                        if isequal(measure,'diceFromOthers') || isequal(measure,'dicePairs')
                            currentMasks = imgCell{tt}{rr}.masks;
                            currentMasks(k) = [];
                            %removes annotators with 0 seniority from the mean
                            masksForMeanCalculation = imgCell{tt}{rr}.masks;
                            indToRemove = [k];
                            for ll=1:length(imgCell{tt}{rr}.masks)
                                [~,currName,~] = fileparts(imgCell{tt}{rr}.masksFilesNames{ll});
                                if SUREExperimentCls.getRadiologistSeniority(SUREExperimentCls.getRadiologistName(currName)) == 1
                                    indToRemove = [indToRemove,ll];
                                end
                            end
                            masksForMeanCalculation(indToRemove) = [];
                            meanShape = stapleWrapper(currentMasks) > 0.5;
                            %meanShape = Utils.calcSegCellSum(masksForMeanCalculation) > length(masksForMeanCalculation)/2;
                            if isequal(measure,'diceFromOthers')
                                meanVal = meanVal + Unitest.dice(meanShape,imgCell{tt}{rr}.masks{k});
                            elseif isequal(measure,'jaccardPairs')
                                meanVal = meanVal + VariabilityExperiment.getAvgAgreement(imgCell{tt}{rr}.masks{k},masksForMeanCalculation);
                            end
                        end
                    end
                    %name = imgCell{tt}{rr}.imgFileName;
                    if exist('groupNames','var')
                        name = [groupNames{tt}(1:2),num2str(rr)];
                    else
                        name = [num2str(rr),num2str(rr)];
                    end
                    
                    meanVal = meanVal / length(imgCell{tt}{rr}.masks);
                    mapObj(name) = meanVal;
                    groupIndsMap(name) = tt;
                end
            end
            mapAvg = mapObj;
            VariabilityExperiment.displayMapAsBar(mapAvg, groupIndsMap, true, true,groupNames);
            title('average dice result with mean annotation per case');
            xlabel('caes number');
            ylabel('average dice result');
        end
        
        function p = pearson(x,y)
            C = cov(x,y);
            p = C(2)/(std(x)*std(y)); 
        end
        
        function res = calcAverageMapPearsonCoeffieicnt(mapsStrct)
            nMaps = length(mapsStrct);
            labels = keys(mapsStrct{1});
            nAnnotators = length(labels);
            X = zeros(nMaps,nAnnotators);
            for l=1:nAnnotators
                curLabel = labels{l};
                for z=1:nMaps
                    X(z,l)= mapsStrct{z}(curLabel);
                end
            end
            pearsonVals = [];
            for ii=1:size(X,1)
                rii = X(ii,:);
                for jj=ii+1:size(X,1)
                   rjj = X(jj,:);
                   p = VariabilityExperiment.pearson(rii,rjj);
                   pearsonVals = [pearsonVals(:); p];
                end
            end
            res = mean(pearsonVals);
        end
        
        function res = calcAverageNumberOfMapShifts(mapsStrct)
            
            nMaps = length(mapsStrct);
            labels = keys(mapsStrct{1});
            nAnnotators = length(labels);
            
            sumOfDiffs = 0;
            for l=1:nAnnotators
                curRanks = zeros(length(mapsStrct),1);
                curLabel = labels{l};
                for z=1:nMaps
                    sortedMapVals = sort(cell2mat(mapsStrct{z}.values));
                    curRanks(z) = find(mapsStrct{z}(curLabel)==sortedMapVals,1);
                end
                curMean = mean(curRanks);
                sumOfDiffs = sumOfDiffs + sum(abs(curMean-curRanks))/nMaps;
            end
            res = sumOfDiffs / nAnnotators;
        end
        
        function res = calcMapMean(mapsStrct)
            
            nMaps = length(mapsStrct);
            labels = keys(mapsStrct{1});
            nAnnotators = length(labels);
            
            res = containers.Map();
            for l=1:nAnnotators
                curRanks = zeros(length(mapsStrct),1);
                curLabel = labels{l};
                for z=1:nMaps
                    sortedMapVals = sort(cell2mat(mapsStrct{z}.values));
                    curRanks(z) = find(mapsStrct{z}(curLabel)==sortedMapVals,1);
                end
                curMean = mean(curRanks);
                res(curLabel) = curMean;
            end
        end
        
        
        function [mapAvg,mapCell] = rankRadiologistsVar(imgCell, measure, sortBySeniority)
            %[possible,consencus,variaiblity] = Utils.calcUnionIntersection(currentMasks);
            if ~exist('measure','var')
                measure = 'diceFromOthers';
            end
            if ~exist('sortBySeniority','var')
                sortBySeniority = true;
            end
            [m,n] = VariabilityExperiment.getSubPlotDim(length(imgCell)+1);
            m=3;n=3;
            mapObj = containers.Map();
            mapCell = {};
            for tt=1:length(imgCell)
                fprintf('case name jacFromMean addToVar\n');
                [~,caseName,~] = fileparts(imgCell{tt}.imgFileName);
                %[~,~,totalVar] = Utils.calcUnionIntersection(imgCell{tt}.masks);
                varChangeRatVec = zeros(length(imgCell{tt}.masks),1);
                diceFromMeanVec =  zeros(length(imgCell{tt}.masks),1);
                for k=1:length(imgCell{tt}.masks)
                    [~,name,~] = fileparts(imgCell{tt}.masksFilesNames{k});
                    [name] = SUREExperimentCls.getRadiologistName(name);
                    name
                    if isequal(measure,'diceFromOthers') || isequal(measure,'dicePairs')
                        currentMasks = imgCell{tt}.masks;
                        currentMasks(k) = [];
                        %removes annotators with 0 seniority from the mean
                        masksForMeanCalculation = imgCell{tt}.masks;
                        indToRemove = [k];
                        for ll=1:length(imgCell{tt}.masks)
                            [~,currName,~] = fileparts(imgCell{tt}.masksFilesNames{ll});
                            if SUREExperimentCls.getRadiologistSeniority(SUREExperimentCls.getRadiologistName(currName)) == 1
                                indToRemove = [indToRemove,ll];
                                SUREExperimentCls.getRadiologistName(currName)
                            end
                        end
                        masksForMeanCalculation(indToRemove) = [];
                        meanShape = stapleWrapper(currentMasks) > 0.5;
                        %meanShape = Utils.calcSegCellSum(masksForMeanCalculation) > length(masksForMeanCalculation)/2;
                        if isequal(measure,'diceFromOthers')
                            val = Unitest.dice(meanShape,imgCell{tt}.masks{k});
                        elseif isequal(measure,'jaccardPairs')
                            val = VariabilityExperiment.getAvgAgreement(imgCell{tt}.masks{k},masksForMeanCalculation);
                        end
                        
                        
                    elseif isequal(measure,'size')
                        val = sum(imgCell{tt}.masks{k}(:));
                    else
                        error 'no param chosen'
                    end
                    
                                     
                    mapObj(name) = val;
                    %fprintf('%s %s %f %f\n',caseName,name,jacFromMeanVec(k),varChangeRatVec(k));
                    
                end
                mapCell{tt} = containers.Map(keys(mapObj),values(mapObj));
                subplot(m,n,tt)
 
                VariabilityExperiment.displayMapAsBar(mapObj,VariabilityExperiment.getSeniorityMap(keys(mapObj)),sortBySeniority,sortBySeniority);
                %ylim([0.6,1])
                fprintf('%s . MeanJacFromMean-%f MeanAddToVar-%f STDjacFromMean-%f  STDaddToVar-%f\n',caseName,mean(diceFromMeanVec),mean(varChangeRatVec),std(diceFromMeanVec),std(varChangeRatVec));
                title(caseName)
            end
            subplot(m,n,tt+1)
            mapAvg = VariabilityExperiment.getMapAvg(imgCell,mapCell);
            VariabilityExperiment.displayMapAsBar(mapAvg,VariabilityExperiment.getSeniorityMap(keys(mapAvg)),sortBySeniority,sortBySeniority);
           
            %VariabilityExperiment.displayMapAsBar(mapAvg, sortBySeniority,sortBySeniority);
            %ylim([0.6,1])
            title('mean');
            annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, ...
                'String','dice result with mean annotation per Radiologist', 'EdgeColor', 'none','HorizontalAlignment', 'center');
        end
        
        function meanAreas = getMeanAreaForImcell(imCell,notScalePixDim)
            meanAreas = zeros(length(imCell),1);
            
            if exist('notScalePixDim','var') && ~notScalePixDim
                hasDims = false;
            else
                hasDims = ~isempty(imCell{1}.dimensions);
            end
            
            for t=1:length(imCell)
                if hasDims
                    multFactor =  imCell{t}.dimensions(1)*imCell{t}.dimensions(2)*imCell{t}.dimensions(3);
                else
                    multFactor = 1;
                end
                curSegAreas = zeros(1,length(imCell{t}.masks));
                for s=1:length(imCell{t}.masks)
                    curSegAreas(s) = sum(imCell{t}.masks{s}(:))*multFactor;
                end
                meanAreas(t) = mean(curSegAreas);
            end
            
        end
        
        %conPosSD - consensus possible surface distance
        function [variability,consensus,possible,allBinVec, conPosSD] = calcObserverVariabilityMeasures(masks,nObservers,dims,calcSD)
            %generates all possible combinations
            nSegs = length(masks);
            allBinVec = logical(dec2bin(1:2^(nSegs)-1)-48);
            allBinVec = logical(allBinVec(sum(allBinVec,2)==nObservers,:));
            
            for z=1:size(allBinVec,1)
                currentMasks = masks(allBinVec(z,:));
                [curPossible,curConsensus,curVariability] = Utils.calcUnionIntersection(currentMasks);
                if  exist('dims','var') && ~isempty(dims)
                    multFactor = dims(1)*dims(2)*dims(3);
                    multFactor2D = dims(1)*dims(2);
                end
                variability(z) = sum(curVariability(:))*multFactor;
                consensus(z) = sum(curConsensus(:))*multFactor;
                possible(z) = sum(curPossible(:))*multFactor;
                
                %multFactor is caluclated inside the function
                if exist('calcSD','var') && calcSD
                    conPosSD(z)  = Utils.surfaceDistanceWithDim(curPossible,curConsensus,dims);
                end
            end
            
            
        end
        
        function [variability2,variability3,variability4] = calcSecondOrderVariabilityMeasures(masks,allBinVec,dims)
            
            for z=1:size(allBinVec,1)
                currentMasks = masks(allBinVec(z,:));
                curSum = Utils.calcSegCellSum(currentMasks);
                N = length(currentMasks);
                curVar2 = curSum >=2 & curSum <= N-2;
                curVar3 = curSum >=3 & curSum <= N-3;
                curVar4 = curSum >=4 & curSum <= N-4;
                if  exist('dims','var') && ~isempty(dims)
                    multFactor = dims(1)*dims(2)*dims(3);
                end
                variability2(z) = sum(curVar2(:))*multFactor;
                variability3(z) = sum(curVar3(:))*multFactor;
                variability4(z) = sum(curVar4(:))*multFactor;
            end
        end
        
        function generateImgStrctsCSVFile(imgStrcts)
            %str = '';
            for z=1:length(imgStrcts)
                VariabilityExperiment.generateImgStrctCSVFile(imgStrcts{z});
                %    currentStr = VariabilityExperiment.generateImgStrctCSVstr(imgStrcts{z});
                %    str = [str, currentStr];
            end
            %fileID = fopen(fileName,'w');
            %fprintf(fileID,str);
            %fclose(fileID);
            
        end
        
        function [meanOut, fullOut] = getCSVFilesRes(csvFilesPath, ignoreNames)
            if ~exist('csvFilesPath','var')
                [csvFilesNames,csvFilesDir] = uigetfile('multiselect','on',{'SUREGrandExp\csvRes\*.csv'});
                if isstr(csvFilesNames)
                    csvFilesNames = {csvFilesNames};
                end
                for z=1:length(csvFilesNames)
                    csvFilesPath{z} = [csvFilesDir,csvFilesNames{z}];
                end
            end
            if ~exist('ignoreNames','var')
                ignoreNames = {};
            end
            if ~iscell(csvFilesPath)
                error 'expected for cell input';
            end
            for z=1:length(csvFilesPath)
                fid = fopen(csvFilesPath{z});
                tline = fgetl(fid);
                header = strsplit(tline,',');
                fclose(fid);
                if ~isequal(header{1},'variability') || ~isequal(header{2},'possible')  || ~isequal(header{3},'consensus') ...
                        || ~isequal(header{4},'meanSegVol')
                    error 'invalid csv file';
                end
                
                csvRes = csvread(csvFilesPath{z},1,0);
                indsToIgnore = [];
                for hh=1:length(header)
                    for kk=1:length(ignoreNames)
                        if ~isempty(findstr(header{hh},ignoreNames{kk}))
                            indsToIgnore = [indsToIgnore,kk];
                            fprintf 'ignoring: header{z} \n';
                        end
                    end
                end
                csvRes(:,indsToIgnore) = [];
                
                radiologistsIndicator = csvRes(:,4:end);
                variability = csvRes(:,1);
                possible = csvRes(:,2);
                consensus = csvRes(:,3);
                nRadiologists = size(radiologistsIndicator,2);
                for t=1:nRadiologists
                    currentRows = sum(radiologistsIndicator,2)==t;
                    fullOut.varVol(t,z) = mean(variability(currentRows,:));
                    fullOut.consVol(t,z) =  mean(consensus(currentRows,:));
                    fullOut.posVol(t,z) = mean(possible(currentRows,:));
                    fullOut.varVolMin(t,z) = min(variability(currentRows,:));
                    fullOut.consVolMin(t,z) =  min(consensus(currentRows,:));
                    fullOut.posVolMin(t,z) = min(possible(currentRows,:));
                    fullOut.varVolMax(t,z) = max(variability(currentRows,:));
                    fullOut.consVolMax(t,z) =  max(consensus(currentRows,:));
                    fullOut.posVolMax(t,z) = max(possible(currentRows,:));
                end
                
            end
            meanOut.varVol = mean(fullOut.varVol,2);
            meanOut.consVol = mean(fullOut.consVol,2);
            meanOut.posVol = mean(fullOut.posVol,2);
            meanOut.varVolMin = mean(fullOut.varVolMin,2);
            meanOut.consVolMin = mean(fullOut.consVolMin,2);
            meanOut.posVolMin = mean(fullOut.posVolMin,2);
            meanOut.varVolMax = mean(fullOut.varVolMax,2);
            meanOut.consVolMax = mean(fullOut.consVolMax,2);
            meanOut.posVolMax = mean(fullOut.posVolMax,2);
        end
        
        function generateImgStrctCSVFile(imgStrct)
            folder = 'SUREGrandExp\csvRes\';
            [~, imageName] =  fileparts(imgStrct.imgFileName);
            frameNum = strrep(num2str(imgStrct.frame),'  ',',');
            fileName = [imageName, '_' frameNum];
            str = VariabilityExperiment.generateImgStrctCSVstr(imgStrct);
            fileID = fopen([folder fileName  '.csv'],'w');
            fprintf(fileID,str);
            fclose(fileID);
        end
        function str = generateImgStrctCSVstr(imgStrct)
            str = ['variability,possible,consensus'];
            for z=1:length(imgStrct.masksFilesNames)
                [~,name] = fileparts(imgStrct.masksFilesNames{z});
                str = [str, ',' name];
            end
            str = [str '\n'];
            [~, imageName] =  fileparts(imgStrct.imgFileName);
            frameNum = strrep(num2str(imgStrct.frame),'  ',',');
            
            nSegs = length(imgStrct.masks);
            for nObservers=1:nSegs
                [variability,consensus,possible,allBinVec] = VariabilityExperiment.calcObserverVariabilityMeasures(imgStrct.masks,nObservers,imgStrct.dimensions);
                for z=1:size(allBinVec,1)
                    currentLine = [num2str(variability(z))];
                    currentLine = [currentLine ',' num2str(possible(z))];
                    currentLine = [currentLine ',' num2str(consensus(z))];
                    currentLine = [currentLine ',' strrep(num2str(allBinVec(z,:)),'  ',',')];
                    str = [str currentLine '\n'];
                end
            end
        end
        
        
        
        function [varVol, varVolMin, varVolMax] = uncertaintyAndCertaintyPlots(imCell,datasetTypeStr)
            if ~exist('datasetTypeStr','var')
                datasetTypeStr = '';
            else
                datasetTypeStr = [': ' datasetTypeStr];
            end
            
            
            N = length(imCell);
            minNSegs = min(inf,length(imCell{1}.masks));
            for t=1:N
                %t
                nSegs = length(imCell{t}.masks);
                minNSegs = min(nSegs,minNSegs);
                for z=1:nSegs
                    %z
                    [variability,consensus,possible] = VariabilityExperiment.calcObserverVariabilityMeasures(imCell{t}.masks,z,imCell{t}.dimensions);
                    varVol(t,z) = mean(variability);
                    consVol(t,z) =  mean(consensus);
                    posVol(t,z) = mean(possible);
                    varVolMin(t,z) = min(variability);
                    consVolMin(t,z) =  min(consensus);
                    posVolMin(t,z) = min(possible);
                    varVolMax(t,z) = max(variability);
                    consVolMax(t,z) =  max(consensus);
                    posVolMax(t,z) = max(possible);
                end
            end
            
            if size(consVol,1)>1
                consVol = mean(consVol);
                varVol = mean(varVol);
                posVol = mean(posVol);
                varVolMin = mean(varVolMin);
                consVolMin = mean(consVolMin);
                posVolMin = mean(posVolMin);
                varVolMax = mean(varVolMax);
                consVolMax = mean(consVolMax);
                posVolMax = mean(posVolMax);
            end
            
            meanAreas = VariabilityExperiment.getMeanAreaForImcell(imCell);
            meanSegAreas = mean(meanAreas);
            
            
            figure, hold on;
            set(gca,'XTick',1:max(size((consVol))) );
            plot(consVol,'-o');
            plot(posVol,'-o');
            plot(ones(size(consVol))*meanSegAreas,'-.');
            line([1:minNSegs;1:minNSegs],[consVolMin;consVolMax],'Color','black')
            line([1:minNSegs;1:minNSegs],[posVolMin;posVolMax],'Color','black')
            xlabel('num of annotators');
            ylabel('area');
            legend('MIN','MAX','mean-seg-area(N)');
            title(['MIN and MAX area as function of number of annotators' datasetTypeStr]);
            hold off;
            
            figure, hold on;
            set(gca,'XTick',1:max(size((consVol))) );
            plot(varVol,'-o');
            line([1:minNSegs;1:minNSegs],[varVolMin;varVolMax],'Color','black')
            xlabel('num of annutators');
            ylabel('area');
            title(['VAR area as function of number of annotators' datasetTypeStr]);
            
            
            figure, hold on;
            set(gca,'XTick',1:max(size((consVol))) );
            plot(varVol/meanSegAreas,'-o');
            line([1:minNSegs;1:minNSegs],[varVolMin/meanSegAreas;varVolMax/meanSegAreas],'Color','black')
            xlabel('num of annutators');
            ylabel('area');
            title(['VAR-seg-area/mean-seg-area as function of number of annotators' datasetTypeStr]);
        end
        
        function resStrct = holdUncExperiment(imCell,intensPrSize,intensityT,outDir)
            
            if exist('outDir','var') && ~exist(outDir,'dir')
                mkdir(outDir);
            end
            
            N = length(imCell);
            resStrct = struct;
            for t=1:N
                seg = imCell{t}.masks{1};
                segBound = Utils.getBoundries(seg);
                gt = segBound & imCell{t}.uncMask;
                
                [~,~,uncertaintyMask] = Utils.calcUnionIntersection(imCell{t}.masks);
                %im = mat2gray(imCell{t}.img);
                
                
                im = Utils.optimizeImContrast2(imCell{t}.img,imCell{t}.masks{1});
                
                %defines params
                params.intensityLocal = true;
                params.filterSize = intensPrSize;
                params.min = true;
                %params.minUncertaintyRegionSize = minUncSize;
                params.intensityLocalT = VariabilityExperiment.calculateThreshold(seg,im);
                [pr] = Prior.intensityPriorLocal(im,seg,intensPrSize,intensityT);
                
                res = pr < intensityT;
                
                if exist('outDir','var')
                    close all;
                    LineDisplay.displayVariability(im,seg,pr, imCell{t}.uncMask);
                    im0 = LineDisplay.getCroppedFrameFromFigure();
                    close all;
                    LineDisplay.displaySegsOverlay(im,zeros(size(im)));
                    resizedIm = LineDisplay.getCroppedFrameFromFigure();
                    close all
                    combIm = [resizedIm,im0];
                    
                    outName = [outDir '\' num2str(t) '-' num2str(t) '.png'];
                    imwrite(combIm,outName);
                end
            end
        end
        
        function T = calculateThreshold(seg,im)
            morphSize = 3;
            erodedSeg = imerode(seg,strel('disk',morphSize));

            T = 1.5*std(single(im(erodedSeg)));
        end
        
        function [resStrct,outMasks] = holdExperiment(imCell,intensPrSize,useMeanSegmentation,outDir, winLength,trivialRun, activeContourRun, segNumToUse, ...
                OptionsIn, OptionsOut)
            
            if exist('outDir','var') && ~isempty(outDir) && ~exist('outDir','dir')
                mkdir(outDir);
            end
            if exist('segNumToUse', 'var') && segNumToUse > 0
                useGivenSegmentation = true;
                useMeanSegmentation = false;
            elseif ~exist('useMeanSegmentation','var') || isempty(useMeanSegmentation)
                useGivenSegmentation = false;
                useMeanSegmentation = true;
            end
           
            N = length(imCell);
            outMasks = cell(size(imCell));
            resStrct = struct;
            for t=1:N
%                 t

                if useMeanSegmentation
                    %seg = Utils.calcSegCellSum(imCell{t}.masks) > length(imCell{t}.masks)/2;
                    seg = stapleWrapper(imCell{t}.masks) > 0.5;
                elseif useGivenSegmentation
                    seg = imCell{t}.masks{segNumToUse};
                else
                    seg = VariabilityEstimator.calcMeanShapeMultSegs(imCell{t}.masks);
                end
                [unionMask,intersectionMask,uncertaintyMask] = Utils.calcUnionIntersection(imCell{t}.masks);
                              %calcs T
                im = imCell{t}.img;
                params.intensityLocalT = VariabilityExperiment.calculateThreshold(seg,im);
                %defines params
                params.intensityLocal = true;
                params.filterSize = intensPrSize;
                params.min = true;
                params.minUncertaintyRegionSize = 2;
                if exist('winLength','var')
                    params.winLength = winLength;
                end
                
                
                %evaluates mask
                if exist('trivialRun','var') && trivialRun
                    [pr] = Prior.intensityPriorLocal(im,seg,intensPrSize, params.intensityLocalT);
                    varMaskEstimated = VariabilityEstimator.evaluate3DVarMaskTrivial(im,logical(seg),pr,params,0.5);
                elseif exist('activeContourRun', 'var') && activeContourRun == 1
                    varMaskEstimated = VariabilityEstimator.evaluate3DVarMaskActiveContours(im, seg);
                elseif exist('activeContourRun', 'var') && activeContourRun == 2
                    varMaskEstimated = VariabilityEstimator.evaluate3DVarMaskSnakes(im, seg, OptionsIn, OptionsOut);
                else
                    [pr] = Prior.intensityPriorLocal(im,seg,intensPrSize, params.intensityLocalT);
                    varMaskEstimated = VariabilityEstimator.evaluate3DVarMask(im,logical(seg),pr,params,0.5);
                end
                outMasks{t} = varMaskEstimated;
                segOrVarMask = varMaskEstimated | logical(seg);
                segAndNotVarMask = ~varMaskEstimated & (logical(seg));
                
                %calcs error
                resStrct.var_dice(t) = (2*sum(sum(sum(varMaskEstimated&uncertaintyMask))))/(sum(varMaskEstimated(:))+sum(uncertaintyMask(:)));
                resStrct.max_dice(t) = (2*sum(sum(sum(segOrVarMask&unionMask))))/(sum(segOrVarMask(:))+sum(unionMask(:)));
                resStrct.min_dice(t) = (2*sum(sum(sum(segAndNotVarMask&intersectionMask))))/(sum(segAndNotVarMask(:))+sum(intersectionMask(:)));
                
                
                dim = imCell{t}.dimensions;
                resStrct.var_volGT(t) = sum(uncertaintyMask(:))*prod(dim);
                resStrct.var_vol(t) = sum(varMaskEstimated(:))*prod(dim);
                resStrct.meanPixAreas = VariabilityExperiment.getMeanAreaForImcell(imCell,false);
                
                
                
                innerVolGT = sum(intersectionMask(:))*prod(dim);
                outerVolGT = sum(unionMask(:))*prod(dim);
                volGT = sum(seg(:));
                resStrct.var_percentage_range_gt(:,t) = [(volGT-innerVolGT)/volGT,(outerVolGT-volGT)/volGT];
                resStrct.var_range_gt(:,t) = [innerVolGT,outerVolGT];
                
                innerVol = sum(sum(sum(seg & ~varMaskEstimated)))*prod(dim);
                outerVol = sum(sum(sum(seg | varMaskEstimated)))*prod(dim);
                resStrct.var_percentage_range(:,t) = [(volGT-innerVol)/volGT,(outerVol-volGT)/volGT];
                resStrct.var_range(:,t) = [innerVol,outerVol];
                
                
                recistOrig = Utils.calcRecistWithDim(seg,dim);
                recistPossibleGT = Utils.calcRecistWithDim(unionMask,dim);
                recistConsensusGT =  Utils.calcRecistWithDim(intersectionMask,dim);
                recistPossibleEst = Utils.calcRecistWithDim(segOrVarMask,dim);
                recistConsensusEst =  Utils.calcRecistWithDim(segAndNotVarMask,dim);
                resStrct.recist_percentage_range_gt(:,t) =  [(recistOrig-recistConsensusGT)/recistOrig,(recistPossibleGT-recistOrig)/recistOrig];
                resStrct.recist_range_gt(:,t) =  [recistConsensusGT,recistPossibleGT];
                
                resStrct.recist_percentage_range(:,t) =  [(recistOrig-recistConsensusEst)/recistOrig,(recistPossibleEst-recistOrig)/recistOrig];
                resStrct.recist_range(:,t) =  [recistConsensusEst,recistPossibleEst];
                
                %display variability
                if exist('outDir','var') && ~isempty(outDir)
                    close all;
                    % display only segmentation without uncertainty region
                    LineDisplay.displayVariabilityFromMask(im, seg, varMaskEstimated, false);
                    im0 = LineDisplay.getCroppedFrameFromFigure();
                    close all;
                    % display segmentation and estimated uncertainty region
                    LineDisplay.displayVariabilityFromMask(im, seg, varMaskEstimated, true);
                    im1 = LineDisplay.getCroppedFrameFromFigure();
                    close all;
                    % display segmentation and ground truth uncertainty region
                    LineDisplay.displaySegsOverlay(im,seg,uncertaintyMask);
                    im2 = LineDisplay.getCroppedFrameFromFigure();
                    close all
                    % display only image without segmentation or uncertainty region
                    LineDisplay.displaySegsOverlay(im,zeros(size(im))); 
                    resizedIm = LineDisplay.getCroppedFrameFromFigure();
                    close all
                    combIm = [resizedIm,im0,im1,im2];
                    outName = [outDir '\' num2str(t) '-' num2str(resStrct.var_dice(t)) '.png'];
                    imwrite(combIm,outName);
                    
                    close all;
                    LineDisplay.displaySegsOverlay(im,seg);
                    imSeg = LineDisplay.getCroppedFrameFromFigure();
                    outSegName = [outDir '\' num2str(t) '-seg' '.png'];
                    outTiffNameOrig = [outDir '\' num2str(t) '-origVar' '.tiff'];
                    outTiffNameEval= [outDir '\' num2str(t) '-evalVar' '.tiff'];
                    outImNameEval= [outDir '\' num2str(t) '-im' '.tiff'];
                    
                    LineDisplay.saveVariabilityAsTiff(im, seg, [],uncertaintyMask,outTiffNameOrig);
                    if exist('activeContourRun', 'var') && activeContourRun > 0
                        LineDisplay.saveVariabilityAsTiff(im, seg, [],varMaskEstimated,outTiffNameEval);
                    else
                        LineDisplay.saveVariabilityAsTiff(im, seg, pr,varMaskEstimated,outTiffNameEval);
                    end
                    IO.saveAsTiff(outImNameEval,im);
                    
                    imwrite(imSeg,outSegName);
                    
                    close all;
                    
                end
            end
            
        end
    end
    
end

