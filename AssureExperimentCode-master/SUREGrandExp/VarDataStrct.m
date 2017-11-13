classdef VarDataStrct
    %VARDATASTRCT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        surfaceDistanceConPos;
        variability;
        possible
        consensus,
        radiologistsNames,
        radiologists,
        imageName,
        frameNum,
        meanSegSize,
        meanRECIST
        
    end
    
    methods (Static)
        
        function [minRes,maxRes, meanRes] = getMinMaxMeanFieldForStrcts(varDataStrcts,fieldName)
            
            nAnnotators = inf;
            for t=1:length(varDataStrcts)
                nAnnotators = min(nAnnotators,size(varDataStrcts{t}.radiologists,2));
            end
            minRes = zeros(nAnnotators,length(varDataStrcts));
            maxRes = zeros(nAnnotators,length(varDataStrcts));
            meanRes = zeros(nAnnotators,length(varDataStrcts));
            for t=1:length(varDataStrcts)
                [curMin,curMax,curMean] = varDataStrcts{t}.getMinMaxMeanField(fieldName);
                minRes(:,t) = curMin(1:nAnnotators);
                maxRes(:,t) = curMax(1:nAnnotators);
                meanRes(:,t) = curMean(1:nAnnotators);
            end
            
        end
        
        
        function generatePlots()
            
        end
        
        function plotAll(varDataStrcts)
            figure
            ax1 = subplot(2,2,1);
            VarDataStrct.plotVarDivArea(varDataStrcts)
            
            ax2 = subplot(2,2,2);
            VarDataStrct.plotVarArea(varDataStrcts);
            
            ax3 = subplot(2,2,3);
            VarDataStrct.plotMinMax(varDataStrcts);
            
            ax4 = subplot(2,2,4);
            VarDataStrct.plotVarDivVarStar({varDataStrcts});
        end
        
        function printMinMaxVarForStrcts(varCells, N)
            
            for z=1:length(varCells)
                
                dataStrctCell = varCells{z};
                meanSeg = cellfun(@(x)x.meanSegSize,dataStrctCell);
                fprintf('caseNum '); fprintf('%i ', 1:length(dataStrctCell)'); fprintf('\n');
                [consVolMin,consVolMax, consVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(dataStrctCell,'consensus');
                [posVolMin,posVolMax, posVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(dataStrctCell,'possible');
                [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(dataStrctCell,'variability');
                
                consVolMin = round(consVolMin ./ repmat(meanSeg',size(consVolMin,1),1),3);
                consVolMax = round(consVolMax ./ repmat(meanSeg',size(consVolMax,1),1),3);
                consVol = round(consVol ./ repmat(meanSeg',size(consVol,1),1),3);
                posVolMin = round(posVolMin ./ repmat(meanSeg',size(posVolMin,1),1),3);
                posVolMax = round(posVolMax ./ repmat(meanSeg',size(posVolMax,1),1),3);
                posVol = round(posVol ./ repmat(meanSeg',size(posVol,1),1),3);
                varVolMin = round(varVolMin ./ repmat(meanSeg',size(varVolMin,1),1),3);
                varVolMax = round(varVolMax ./ repmat(meanSeg',size(varVolMax,1),1),3);
                varVol = round(varVol ./ repmat(meanSeg',size(varVol,1),1),3);
                if ~exist('N','var')
                    fprintf('%s ',consVolMin);  fprintf('\n');
                    fprintf('%s ', consVolMax); fprintf('\n');
                    fprintf('%s ', consVol); fprintf('\n');
                    
                    fprintf('%s ', posVolMin); fprintf('\n');
                    fprintf('%s ', posVolMax); fprintf('\n');
                    fprintf('%s ', posVol); fprintf('\n');
                    
                    fprintf('%s ', varVolMin); fprintf('\n');
                    fprintf('%s ', varVolMax); fprintf('\n');
                    fprintf('%s ', varVol); fprintf('\n');
                else
                    if N(z)<size(consVolMin,1)
                        fprintf('consensus(min) '); fprintf('%f ',consVolMin(N(z),:));  fprintf('\n');
                        fprintf('consensus(max) '); fprintf('%f ', consVolMax(N(z),:));  fprintf('\n');
                    end
                    
                    fprintf('consensus '); fprintf('%f ', consVol(N(z),:));  fprintf('\n');
                    
                    if N(z)<size(consVolMin,1)
                        fprintf('possible(min) '); fprintf('%f ', posVolMin(N(z),:));   fprintf('\n');
                        fprintf('possible(max) '); fprintf('%f ', posVolMax(N(z),:)); fprintf('\n');
                    end
                    fprintf('possible'); fprintf(' %f ', posVol(N(z),:));  fprintf('\n');
                    
                    if N(z)<size(consVolMin,1)
                        fprintf('var(min)'); fprintf(' %f ', varVolMin(N(z),:)); fprintf('\n');
                        fprintf('var(max)'); fprintf(' %f ', varVolMax(N(z),:)); fprintf('\n');
                    end
                    fprintf('var'); fprintf(' %f ', varVol(N(z),:)); fprintf('\n');
                    fprintf('\n\n');
                end
            end
            
        end
        
        function plotVarDivArea(varDataStrcts)
            [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability');
            varVolMaxDiv = varVolMax;
            varVolMinDiv = varVolMin;
            varVolDiv = varVol;
            for t=1:length(varDataStrcts)
                varVolMaxDiv(:,t) = varVolMax(:,t) / varDataStrcts{t}.meanSegSize;
                varVolMinDiv(:,t) = varVolMin(:,t) / varDataStrcts{t}.meanSegSize;
                varVolDiv(:,t) = varVolDiv(:,t) / varDataStrcts{t}.meanSegSize;
            end
            varVolMaxDiv = mean(varVolMaxDiv,2);
            varVolMinDiv = mean(varVolMinDiv,2);
            varVolDiv =  mean(varVolDiv,2);
            nRadiologists = size(varVolMax,1);
            
            hold on;
            set(gca,'XTick',1:nRadiologists );
            plot(varVolDiv,'-o');
            for z=1:length(varVolDiv)
                text(z+0.1,varVolDiv(z)+0.005,num2str(varVolDiv(z),2));
            end
            line([(1:nRadiologists);(1:nRadiologists)],[varVolMinDiv';varVolMaxDiv'],'Color','black')
            xlabel('num of annutators');
            ylabel('VAR-seg-area/mean-seg-area');
            title(['VAR-seg-area/mean-seg-area as function of number of annotators ']);
            hold off;
        end
        
        function plotVarArea(varDataStrcts, plotSecOrder)
            if ~exist('plotSecOrder')
                plotSecOrder = false;
            end
            [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability');
            varVolMin = mean(varVolMin,2);
            varVolMax = mean(varVolMax,2);
            varVol = mean(varVol,2);
            nRadiologists = size(varVolMax,1);
            
            
            %var
            hold on;
            set(gca,'XTick',1:nRadiologists);
            plot(varVol,'-o');
            line([1:nRadiologists;1:nRadiologists],[varVolMin';varVolMax'],'Color','black')
            xlabel('num of annutators');
            ylabel('area');
            title(['VAR area as function of number of annotators ']);
            
            if plotSecOrder
                [var2VolMin,var2VolMax, varVol2] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability2');
                
                var2VolMin = mean(var2VolMin,2);
                var2VolMax = mean(var2VolMax,2);
                varVol2 = mean(varVol2,2);
                plot(varVol2,'-o','Color','red');
                line([1:nRadiologists;1:nRadiologists],[var2VolMin';var2VolMax'],'Color','red');
                
                [var3VolMin,var3VolMax, varVol3] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability3');
                var3VolMin = mean(var3VolMin,2);
                var3VolMax = mean(var3VolMax,2);
                varVol3 = mean(varVol3,2);
                plot(varVol3,'-o','Color','blue');
                line([1:nRadiologists;1:nRadiologists],[var3VolMin';var3VolMax'],'Color','blue');
                
                [var4VolMin,var4VolMax, varVol4] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability3');
                var4VolMin = mean(var4VolMin,2);
                var4VolMax = mean(var4VolMax,2);
                varVol4 = mean(varVol4,2);
                plot(varVol4,'-o','Color','blue');
                line([1:nRadiologists;1:nRadiologists],[var4VolMin';var4VolMax'],'Color','blue');
            end
            hold off;
            
        end
        
        
        
        
        
        function [varVolCell,varVolMinCell,varVolMaxCell]= plotVarDivVarStar(varDataStrcts, titles, divByStar, calcSD)
            if ~exist('calcSD','var')
                calcSD = false;
            end
            
            if ~iscell(varDataStrcts{1})
                error 'function receives cell of cells'
            end
            
            if ~exist('divByStar','var')
                divByStar = true;
            end
            varVolMinCell = {};
            varVolMaxCell = {};
            varVolCell = {};
            
            nRadiologists = inf;
            namesCell = {};
            for t=1:length(varDataStrcts)
                [~,~, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability');
                nRadiologists = min(nRadiologists,size(varVol,1));
                if exist('titles','var') && ~isempty(titles)
                    nextName = titles{t};
                else
                    curNames = '';
                    for k=1:length(varDataStrcts{t})
                        curNames = [curNames ', ' varDataStrcts{t}{k}.imageName];
                    end
                    nextName = curNames;
                end
                namesCell = {namesCell{:},[strrep(nextName,'_','-')]};
            end
            
            
            for t=1:length(varDataStrcts)
                if ~calcSD
                    [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability');
                    for k=1:size(varVol,2)
                        if divByStar
                            varVol(:,k) = varVol(:,k) / varVol(nRadiologists,k);
                            varVolMin(:,k) = varVolMin(:,k) /varVolMin(nRadiologists,k);
                            varVolMax(:,k) = varVolMax(:,k) / varVolMax(nRadiologists,k);
                            N = nRadiologists;
                        else
                            varVol(:,k) = varVol(:,k) / varDataStrcts{t}{k}.meanSegSize;
                            varVolMin(:,k) = varVolMin(:,k) / varDataStrcts{t}{k}.meanSegSize;
                            varVolMax(:,k) = varVolMax(:,k) / varDataStrcts{t}{k}.meanSegSize;
                            N = size(varVolMax,1);
                        end
                    end
                else
                    [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'surfaceDistanceConPos');
                        
                    for k=1:size(varVol,2)
                    %todo -ugly, should use correct names in the surface
                    %distance case.. but didnt have time
                        varVol(:,k) = varVol(:,k) / varVol(nRadiologists,k);
                        varVolMin(:,k) = varVolMin(:,k) /varVolMin(nRadiologists,k);
                        varVolMax(:,k) = varVolMax(:,k) / varVolMax(nRadiologists,k);
                        N = nRadiologists;
                    end
                end
                varVolMinCell{t} = mean(varVolMin(1:N,:),2);
                varVolMaxCell{t} = mean(varVolMax(1:N,:),2);
                varVolCell{t} = mean(varVol(1:N,:),2);
                
                
            end
            
            
            %var
            hold on;
            colors = Utils.getRGBColors();
            colors(1:4,:) = [44,118,167 ; 105,141,56;210,87,34;92,196,239]/255;
            set(gca,'XTick',1:nRadiologists);
            for rr=1:length(varDataStrcts)
                plot(varVolCell{rr},'-o','color',colors(rr,:));
                
                
            end
            if length(varDataStrcts)==1
                for z=1:length(varVolCell{1})
                    text(z+0.1,varVolCell{1}(z)+0.005,num2str(varVolCell{1}(z),2));
                end
            end
            legend(namesCell);
            %legend('liver tumors','lung tumors','kidneys','liver tumors(min and max bounds)','lung tumors(min and max bounds)','kidneys(min and max bounds)');
            
            for rr=1:length(varDataStrcts)
                plot(varVolMinCell{rr},'.','color',colors(rr,:),'MarkerSize', 12);
                plot(varVolMaxCell{rr},'.','color',colors(rr,:),'MarkerSize', 12);
                
                plotArea = length(varDataStrcts)<5 && (~divByStar);
                if plotArea
                    x1 = 1:length(varVolMinCell{rr});
                    y1 = varVolMinCell{rr}';
                    y2 = varVolCell{rr}';
                    X=[x1,fliplr(x1)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    h = fill(X,Y,colors(rr,:));
                    set(h,'facealpha',.1)
                    set(h,'EdgeColor','None');
                    
                    y1 = varVolMaxCell{rr}';             %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    h = fill(X,Y,colors(rr,:));
                    set(h,'facealpha',.1)
                    set(h,'EdgeColor','None');
                end
            end
            
            
            
            xlabel('num of annotators');
            % if (~divByStar)
            %     ylabel('(VAR)/mean(slice Vol) as a function of number of annotators');
            %     title(['(VAR)/mean(slice Vol) as a function of number of annotators']);
            % else
            %     ylabel('(VAR)/(VAR*) as a function of number of annotators');
            %     title(['(VAR)/(VAR*) as a function of number of annotators']);
            % end
            
            hold off;
            if length(titles)==0
                titles = {''};
            end
            
            for t=1:length(titles)
                fprintf('%s ',strrep(titles{t},' ','_'));
                for z=1:length(varVolCell{t})
                    a = num2str(num2str(varVolCell{t}(z),'%.2f'));
                    b = num2str(num2str(varVolMaxCell{t}(z)-varVolCell{t}(z),'%.2f'));
                    c = num2str(num2str(varVolCell{t}(z)-varVolMinCell{t}(z),'%.2f'));
                    fprintf('%s+%s-%s ',a,b,c);
                end
                fprintf('\n',titles{t});
            end
            
        end
        
        function plotSecondOrderVar(varDataStrcts, titles)
            if ~iscell(varDataStrcts{1})
                error 'function receives cell of cells'
            end
            
            
            var2VolMinCell = {}; var3VolMinCell = {}; var4VolMinCell = {};
            var2VolMaxCell = {}; var3VolMaxCell = {}; var4VolMaxCell = {};
            varVol2Cell = {}; varVol3Cell = {}; varVol4Cell = {};
            
            nRadiologists = inf;
            namesCell = {};
            for t=1:length(varDataStrcts)
                [~,~, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability');
                nRadiologists = min(nRadiologists,size(varVol,1));
                
                if exist('titles','var') && ~isempty(titles)
                    nextName = titles{t};
                else
                    curNames = '';
                    for k=1:length(varDataStrcts{t})
                        curNames = [curNames ', ' varDataStrcts{t}{k}.imageName];
                    end
                    nextName = curNames;
                end
                
                %namesCell{t} = strrep(curNames,'_','-');
                
                namesCell = {namesCell{:},[strrep(nextName,'_','-') '-VAR2']};
                namesCell = {namesCell{:},[strrep(nextName,'_','-') '-VAR3']};
                namesCell = {namesCell{:},[strrep(nextName,'_','-') '-VAR4']};
            end
            
            
            for t=1:length(varDataStrcts)
                
                [var2VolMin,var2VolMax, varVol2] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability2DivVariability');
                
                var2VolMin = mean(var2VolMin(1:nRadiologists,:),2);
                var2VolMax = mean(var2VolMax(1:nRadiologists,:),2);
                varVol2 = mean(varVol2(1:nRadiologists,:),2);
                
                [var3VolMin,var3VolMax, varVol3] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability3DivVariability');
                var3VolMin = mean(var3VolMin(1:nRadiologists,:),2);
                var3VolMax = mean(var3VolMax(1:nRadiologists,:),2);
                varVol3 = mean(varVol3(1:nRadiologists,:),2);
                
                [var4VolMin,var4VolMax, varVol4] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts{t},'variability4DivVariability');
                var4VolMin = mean(var4VolMin(1:nRadiologists,:),2);
                var4VolMax = mean(var4VolMax(1:nRadiologists,:),2);
                varVol4 = mean(varVol4(1:nRadiologists,:),2);
                
                varVol2Cell{t} = varVol2;
                var2VolMinCell{t} = var2VolMin;
                var2VolMaxCell{t} = var2VolMax;
                
                varVol3Cell{t} = varVol3;
                var3VolMinCell{t} = var3VolMin;
                var3VolMaxCell{t} = var3VolMax;
                
                varVol4Cell{t} = varVol4;
                var4VolMinCell{t} = var4VolMin;
                var4VolMaxCell{t} = var4VolMax;
            end
            
            
            %var
            hold on;
            colors = Utils.getRGBColors();
            
            set(gca,'XTick',1:nRadiologists);
            for rr=1:length(varDataStrcts)
                plot(varVol2Cell{rr},'--o','color',colors(rr,:));
                plot(varVol3Cell{rr},'-.o','color',colors(rr,:));
                plot(varVol4Cell{rr},'-.o','color',colors(rr,:));
            end
            if length(varDataStrcts)==1
                for z=1:length(varVol2Cell{1})
                    text(z+0.1,varVol2Cell{1}(z)+0.005,num2str(varVol2Cell{1}(z),2));
                    text(z+0.1,varVol3Cell{1}(z)+0.005,num2str(varVol3Cell{1}(z),2));
                    text(z+0.1,varVol4Cell{1}(z)+0.005,num2str(varVol4Cell{1}(z),2))
                end
            end
            legend(namesCell);
            %legend('liver tumors','lung tumors','kidneys','liver tumors(min and max bounds)','lung tumors(min and max bounds)','kidneys(min and max bounds)');
            
            %for rr=1:length(varDataStrcts)
            %    plot(varVolMinCell{rr},'.','color',colors(rr,:),'MarkerSize', 12);
            %    plot(varVolMaxCell{rr},'.','color',colors(rr,:),'MarkerSize', 12);
            %end
            xlabel('num of annotators');
            ylabel('(VAR_k)/(VAR_1) as a function of number of annotators');
            title(['(VAR_k)/(VAR_1) as a function of number of annotators']);
            hold off;
        end
        
        function plotFieldDivValSeveral(varDataStrcts,titles,fieldName,divFieldName)
            
            if ~exist('divFieldName','var')
                divFieldName =  'meanSegSize';
            end
            addText = false;
            plots = {};
            
            for z=1:length(varDataStrcts)
                plots = {plots{:},VarDataStrct.plotFieldDivVal(varDataStrcts{z},fieldName,divFieldName,addText)};
            end
            
            
            xlabel('num of annotators');
            ylabel('');
            legend([plots{:}],titles);
            
            
        end
        
        
        function res = calcPairwiseVarImgCell(imgCell, fieldName,divByMean)
            N = length(imgCell.masks);
            
            
            if ~exist('divByMean','var') || (exist('divByMean','var') && divByMean)
                meanVec = zeros(N,1);
                for z=1:N
                    if strcmp(fieldName,'variability') || strcmp(fieldName,'volDiff')
                        meanVec(z) = sum(imgCell.masks{z}(:))*imgCell.dimensions(1)*imgCell.dimensions(2)*imgCell.dimensions(3);
                    elseif strcmp(fieldName,'surfaceDistanceConPos')
                        meanVec(z) = Utils.calcRecistWithDim(imgCell.masks{z},imgCell.dimensions);
                    end
                end
                divVal = mean(meanVec);
            else
                divVal = 1;
            end
            
            minInd = -1;
            res = zeros((N*(N-1))/2,1);
            t=1;
            for z=1:N
                for k=(z+1):N
                    [~,~,curVariability] = Utils.calcUnionIntersection({imgCell.masks{z},imgCell.masks{k}});
                    if strcmp(fieldName,'variability')
                        res(t) = sum(curVariability(:))*imgCell.dimensions(1)*imgCell.dimensions(2)*imgCell.dimensions(3);
                    elseif strcmp(fieldName,'surfaceDistanceConPos')
                        tempVal = Utils.surfaceDistanceWithDim(imgCell.masks{z},imgCell.masks{k},imgCell.dimensions);
                        if tempVal > min(res)
                            minInd = [z,k];
                        end
                        res(t) = tempVal;
                    elseif  strcmp(fieldName,'volDiff')
                        v1 = imgCell.masks{k}*imgCell.dimensions(1)*imgCell.dimensions(2)*imgCell.dimensions(3);
                        v2 = imgCell.masks{z}*imgCell.dimensions(1)*imgCell.dimensions(2)*imgCell.dimensions(3);
                        v1 = sum(v1(:));
                        v2 = sum(v2(:));
                        res(t) = abs(v1-v2);
                    end
                    t=t+1;
                end
            end
            
            res = res ./ divVal;
        end
        
        
        function plotPairwiseComparison(varStrctCells, fieldName, divFieldName, groupNames, varImgCell)
            if ~iscell(varStrctCells{1})
                error 'expected for cell of cells';
            end
            
            mapObj = containers.Map();
            groupIndsMap =  containers.Map();
            linesMapObj = containers.Map();
            
            
            for tt=1:length(varStrctCells)
                for rr=1:length(varStrctCells{tt})
                    curStrct = varStrctCells{tt}{rr};
                    if exist('varImgCell','var')
                        res = VarDataStrct.calcPairwiseVarImgCell(varImgCell{tt}{rr}, fieldName,~strcmp(divFieldName,'.'));
                        minVal = min(res);
                        maxVal = max(res);
                        val = mean(res);
                    else
                        if exist('divFieldName','var') && ~strcmp(divFieldName,'.')
                            divVal = getfield(curStrct,divFieldName);
                        else
                            divVal = 1;
                        end
                        evaluationVec = getfield(curStrct,fieldName);
                        relevantCases = sum(curStrct.radiologists,2)==2;
                        
                        %meanShape = Utils.calcSegCellSum(masksForMeanCalculation) > length(masksForMeanCalculation)/2;
                        val = mean(evaluationVec(relevantCases))/divVal;
                        minVal = min(evaluationVec(relevantCases))/divVal;
                        maxVal = max(evaluationVec(relevantCases))/divVal;
                    end
                    %name = imgCell{tt}{rr}.imgFileName;
                    if exist('groupNames','var')
                        name = [groupNames{tt}(1:2),num2str(rr)];
                    else
                        name = [num2str(rr),num2str(rr)];
                    end
                    
                    
                    mapObj(name) = val;
                    linesMapObj(name) = [minVal,maxVal];
                    groupIndsMap(name) = tt;
                end
            end
            mapAvg = mapObj;
            VariabilityExperiment.displayMapAsBar(mapAvg, groupIndsMap, true, true,groupNames,linesMapObj);
            title('');
            xlabel('case number');
            ylabel('');
            
        end
        
        function printMinMaxMeanGaps(meanVec, maxVec, minVec, rndStr)
            if ~exist('rndStr','var')
                rndStr = '%.2f';
            end
            for z=1:length(meanVec)
                a = num2str(num2str(meanVec(z),rndStr));
                b = num2str(num2str(maxVec(z)-meanVec(z),rndStr));
                c = num2str(num2str(meanVec(z)-minVec(z),rndStr));
                fprintf('%s+%s-%s ',a,b,c);
            end
            fprintf('\n');
        end
        function [plotRes, minVal,maxVal, meanVal] = plotFieldDivVal(varDataStrcts,fieldName,divFieldName, addText,roundStr)
            if ~exist('roundStr','var')
                roundStr = '%.2f';
            end
            [minVal,maxVal, meanVal] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,fieldName);
            for t=1:length(varDataStrcts)
                if strcmp(divFieldName,'.')
                    divVal = 1;
                else
                    divVal = getfield(varDataStrcts{t},divFieldName);
                end
                
                maxVal(:,t) = maxVal(:,t) / divVal;
                minVal(:,t) = minVal(:,t) / divVal;
                meanVal(:,t) = meanVal(:,t) / divVal;
            end
            
            %maxVal = round(mean(maxVal,2),2);
            %minVal = round(mean(minVal,2),2);
            %meanVal =  round(mean(meanVal,2),2);
            maxVal = mean(maxVal,2);
            minVal = mean(minVal,2);
            meanVal =  mean(meanVal,2);
            nRadiologists = size(meanVal,1);
            
            hold on;
            set(gca,'XTick',1:nRadiologists);
            
            plotRes = plot(meanVal,'-o');
            color = get(plotRes, 'Color');
            %line([1:nRadiologists;1:nRadiologists],[varVolMinDiv';varVolMaxDiv'],'Color','black')
            if ~exist('addText') || (exist('addText') && addText)
                for z=1:length(meanVal)
                    text(z+0.1,meanVal(z)+0.005,num2str(meanVal(z),roundStr));
                end
            end
            
            plot(minVal,'.','color',color,'MarkerSize', 12);
            plot(maxVal,'.','color',color,'MarkerSize', 12);
            
            x1 = 1:length(meanVal);
            y1 = minVal';
            y2 = maxVal';
            X=[x1,fliplr(x1)];                %#create continuous x value array for plotting
            Y=[y1,fliplr(y2)];              %#create y values for out and then back
            h = fill(X,Y,color);
            set(h,'facealpha',.25)
            set(h,'EdgeColor','None');
            
            meanVal(10)
            
        end
        
        
        function plotMinMaxVarDiv(varDataStrcts, titleStr, params)
            if ~exist('titleStr','var')
                titleStr = ['MIN/Mean(slice Vol) and MAX/Mean(slice Vol) as function of number of annotators '];
            end
            if ~exist('params','var')
                params.minMax = true;
                params.var = false;
            end
            
            hold on;
            plots = {};
            if params.var
                plots = {plots{:},VarDataStrct.plotFieldDivVal(varDataStrcts,'variability','meanSegSize')};
            end
            if params.minMax
                plots = {plots{:},VarDataStrct.plotFieldDivVal(varDataStrcts,'consensus','meanSegSize')};
                plots = {plots{:},VarDataStrct.plotFieldDivVal(varDataStrcts,'possible','meanSegSize')};
            end
            
            
            ylimMin = 0.6;  ylimMax = 0.7; legendTitles = {};
            if params.var
                legendTitles = {legendTitles{:},'s-vol(var)/mean annotation s-vol'};
                ylimMin = 0;
            end
            if params.minMax
                legendTitles = {legendTitles{:},'s-vol(consensus)/mean annotation s-vol','s-vol(possible)/mean annotation s-vol'};
                ylimMax = 1.5;
            end
            
            xlabel('num of annotators');
            ylabel('');
            ylim([ylimMin,ylimMax]);
            if length(plots)>1
                legend([plots{:}],legendTitles);
            else
                legend(legendTitles);
            end
            title(titleStr);
            hold off;
            
        end
        
        
        
        function plotMinMax(varDataStrcts)
            [consVolMin,consVolMax, consVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'consensus');
            consVolMin = mean(consVolMin,2);
            consVolMax = mean(consVolMax,2);
            consVol = mean(consVol,2);
            [posVolMin,posVolMax, posVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'possible');
            posVolMin = mean(posVolMin,2);
            posVolMax = mean(posVolMax,2);
            posVol = mean(posVol,2);
            nRadiologists = size(consVolMin,1);
            
            meanSize = 0;
            for z=1:length(varDataStrcts)
                meanSize = meanSize + varDataStrcts{z}.meanSegSize;
            end
            meanSize = meanSize / length(varDataStrcts);
            
            hold on;
            set(gca,'XTick',1:nRadiologists);
            plot(consVol,'-o');
            plot(posVol,'-o');
            plot(ones(size(consVol))*meanSize,'-.');
            line([1:nRadiologists;1:nRadiologists],[consVolMin';consVolMax'],'Color','black')
            line([1:nRadiologists;1:nRadiologists],[posVolMin';posVolMax'],'Color','black')
            xlabel('num of annotators');
            ylabel('area');
            legend('MIN','MAX','mean-seg-area(N)');
            title(['MIN and MAX area as function of number of annotators ']);
            hold off;
        end
        
        function varDataStrctCell = fromImgStrctCell(imgStrctCell)
            varDataStrctCell = cell(length(imgStrctCell),1);
            for z=1:length(imgStrctCell)
                varDataStrctCell{z} = VarDataStrct.fromImgStrct(imgStrctCell{z});
            end
        end
        function varDataStrct = fromImgStrct(imgStrct)
            varDataStrct = VarDataStrct;
            for z=1:length(imgStrct.masksFilesNames)
                [~,name] = fileparts(imgStrct.masksFilesNames{z});
                varDataStrct.radiologistsNames{z} = name;
            end
            [~, varDataStrct.imageName] =  fileparts(imgStrct.imgFileName);
            varDataStrct.frameNum = strrep(num2str(imgStrct.frame),'  ',',');
            
            %meanSegArea
            if ~isempty(imgStrct.dimensions)
                multFactor =  imgStrct.dimensions(1)*imgStrct.dimensions(2)*imgStrct.dimensions(3);
            else
                multFactor = 1;
            end
            
            curSegRECIST = zeros(length(imgStrct.masks),1);
            curSegAreas = zeros(length(imgStrct.masks),1);
            for qq=1:length(imgStrct.masks)
                curSegAreas(qq) = sum(imgStrct.masks{qq}(:))*multFactor;
                
                if ~isempty(imgStrct.dimensions)
                    curSegRECIST(qq) = Utils.calcRecistWithDim(imgStrct.masks{qq},imgStrct.dimensions);
                else
                    curSegRECIST(qq) = Utils.calcRecistWithDim(imgStrct.masks{qq});
                end
            end
            varDataStrct.meanSegSize = mean(curSegAreas);
            varDataStrct.meanRECIST = mean(curSegRECIST);
            nSegs = length(imgStrct.masks);
            varDataStrct.radiologists = [];
            varDataStrct.possible = [];
            varDataStrct.consensus = [];
            
            calcSD = true;
            for nObservers=1:nSegs
                if ~calcSD
                    [var,con,pos,allBinVec] = ...
                        VariabilityExperiment.calcObserverVariabilityMeasures(imgStrct.masks,nObservers,imgStrct.dimensions);
                else
                    [var,con,pos,allBinVec, posConSD] = ...
                        VariabilityExperiment.calcObserverVariabilityMeasures(imgStrct.masks,nObservers,imgStrct.dimensions,true);
                    varDataStrct.surfaceDistanceConPos = [varDataStrct.surfaceDistanceConPos;posConSD'];
                end
                varDataStrct.variability = [varDataStrct.variability;var'];
                varDataStrct.consensus = [varDataStrct.consensus;con'];
                varDataStrct.possible = [varDataStrct.possible;pos'];
                varDataStrct.radiologists = logical([varDataStrct.radiologists;allBinVec]);
                
            end
            
        end
        
        
        
        function plotMinMaxVarDivOld(varDataStrcts, titleStr, params)
            
            if ~exist('titleStr','var')
                titleStr = ['MIN/Mean(slice Vol) and MAX/Mean(slice Vol) as function of number of annotators '];
            end
            
            if ~exist('params','var')
                params.minMax = true;
                params.var = false;
            end
            
            
            if params.minMax
                [consVolMin,consVolMax, consVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'consensus');
                [posVolMin,posVolMax, posVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'possible');
                consVolMaxDiv = consVolMax; posVolMaxDiv = posVolMax;
                consVolMinDiv = consVolMin; posVolMinDiv = posVolMin;
                consVolDiv = consVol; posVolDiv = posVol;
            end
            if params.var
                [varVolMin,varVolMax, varVol] = VarDataStrct.getMinMaxMeanFieldForStrcts(varDataStrcts,'variability');
                varVolMaxDiv = varVolMax;
                varVolMinDiv = varVolMin;
                varVolDiv = varVol;
            end
            
            
            
            for t=1:length(varDataStrcts)
                if params.var
                    varVolMaxDiv(:,t) = varVolMax(:,t) / varDataStrcts{t}.meanSegSize;
                    varVolMinDiv(:,t) = varVolMin(:,t) / varDataStrcts{t}.meanSegSize;
                    varVolDiv(:,t) = varVolDiv(:,t) / varDataStrcts{t}.meanSegSize;
                end
                if params.minMax
                    consVolMaxDiv(:,t) = consVolMax(:,t) / varDataStrcts{t}.meanSegSize;
                    consVolMinDiv(:,t) = consVolMin(:,t) / varDataStrcts{t}.meanSegSize;
                    consVolDiv(:,t) = consVolDiv(:,t) / varDataStrcts{t}.meanSegSize;
                    posVolMaxDiv(:,t) = posVolMax(:,t) / varDataStrcts{t}.meanSegSize;
                    posVolMinDiv(:,t) = posVolMin(:,t) / varDataStrcts{t}.meanSegSize;
                    posVolDiv(:,t) = posVolDiv(:,t) / varDataStrcts{t}.meanSegSize;
                end
            end
            
            if params.minMax
                posVolDiv = mean(posVolDiv,2);
                consVolDiv = mean(consVolDiv,2);
                posVolMinDiv =  mean(posVolMinDiv,2);
                posVolMaxDiv =  mean(posVolMaxDiv,2);
                consVolMinDiv =  mean(consVolMinDiv,2);
                consVolMaxDiv =  mean(consVolMaxDiv,2);
            end
            if params.var
                varVolDiv = mean(varVolDiv,2);
                varVolMaxDiv = mean(varVolMaxDiv,2);
                varVolMinDiv = mean(varVolMinDiv,2);
            end
            nRadiologists = size(consVolMin,1);
            
            
            hold on;
            set(gca,'XTick',1:nRadiologists);
            
            ylimMin = 0.6;  ylimMax = 0.7; legendTitles = {};
            if params.var
                legendTitles = {legendTitles(:),'var/Mean(slice Vol)'};
                ylimMin = 0;
                plot(varVolDiv,'-o');
                %line([1:nRadiologists;1:nRadiologists],[varVolMinDiv';varVolMaxDiv'],'Color','black')
                for z=1:length(varVolDiv)
                    text(z+0.1,varVolDiv(z)+0.005,num2str(varVolDiv(z),'%.2f'));
                end
                
                
                x1 = 1:length(varVolDiv)
                y1 = varVolMinDiv';
                y2 = varVolMaxDiv';
                X=[x1,fliplr(x1)];                %#create continuous x value array for plotting
                Y=[y1,fliplr(y2)];              %#create y values for out and then back
                h = fill(X,Y,'black');
                set(h,'facealpha',.04)
                set(h,'EdgeColor','None');
                
            end
            
            if params.minMax
                legendTitles = {legendTitles{:},'vol(consensus)/vol(annotation)','vol(possible)/vol(annotation)'};
                ylimMax = 1.5;
                plot(consVolDiv,'-o');
                plot(posVolDiv,'-o');
                %line([1:nRadiologists;1:nRadiologists],[consVolMinDiv';consVolMaxDiv'],'Color','black')
                %line([1:nRadiologists;1:nRadiologists],[posVolMinDiv';posVolMaxDiv'],'Color','black')
                for z=1:length(consVolDiv)
                    text(z+0.1,posVolDiv(z)+0.005,num2str(posVolDiv(z),'%.2f'));
                    text(z+0.1,consVolDiv(z)+0.005,num2str(consVolDiv(z),'%.2f'));
                end
                
                x1 = 1:length(consVolDiv)
                y1 = consVolMinDiv';
                y2 = consVolMaxDiv';
                X=[x1,fliplr(x1)];                %#create continuous x value array for plotting
                Y=[y1,fliplr(y2)];              %#create y values for out and then back
                h = fill(X,Y,'black');
                set(h,'facealpha',.04)
                set(h,'EdgeColor','None');
                
                x1 = 1:length(consVolDiv)
                y1 = posVolMinDiv';
                y2 = posVolMaxDiv';
                X=[x1,fliplr(x1)];                %#create continuous x value array for plotting
                Y=[y1,fliplr(y2)];              %#create y values for out and then back
                h = fill(X,Y,'black');
                set(h,'facealpha',.04)
                set(h,'EdgeColor','None');
            end
            
            
            xlabel('num of annotators');
            ylabel('');
            ylim([ylimMin,ylimMax]);
            legend(legendTitles);
            title(titleStr);
            hold off;
            
            
            
        end
    end
    methods
        
        function [minRes,maxRes, meanRes] = getMinMaxMeanField(varDataStrct,fieldName)
            vect = getfield(varDataStrct,fieldName);
            nRadiologists = size(varDataStrct.radiologists,2);
            meanRes = zeros(nRadiologists,1);
            minRes = zeros(nRadiologists,1);
            maxRes = zeros(nRadiologists,1);
            for n=1:nRadiologists
                currentRows = sum(varDataStrct.radiologists,2)==n;
                meanRes(n) = mean(vect(currentRows,:));
                minRes(n) = min(vect(currentRows,:));
                maxRes(n) = max(vect(currentRows,:));
            end
            
        end
    end
    
end

