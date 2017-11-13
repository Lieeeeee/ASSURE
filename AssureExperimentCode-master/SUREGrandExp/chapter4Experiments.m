%annotationsDirs = strcat('C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\',SUREExperimentCls.shortDirNames);
annotationsDirs = strcat('C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\',SUREExperimentCls.novicesDirs);
dataDir = 'C:\Users\drorcohe\Desktop\SUREGrandExperiment\data';
mode3D  = true;

%extracts all liver annotation files 
[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.livers, SUREExperimentCls.windowingLivers);
%generates cell of liver image structs from the short group
liverImgCellNov = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortLivers, mode3D);
%removes novices
liverImgCell = ImgStrct.removeNovices(liverImgCellNov);
%generates data struct cell
liverDataStrctCell = VarDataStrct.fromImgStrctCell(liverImgCell);


[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.lungs, SUREExperimentCls.windowingLungs);
lungsImgCellNov = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortLungs, mode3D);
lungsImgCell = ImgStrct.removeNovices(lungsImgCellNov);
lungDataStrctCell = VarDataStrct.fromImgStrctCell(lungsImgCell);

[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.brains, SUREExperimentCls.windowingBrains);
brainImgCellNov = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortBrains, mode3D);
brainImgCell = ImgStrct.removeNovices(brainImgCellNov);
brainDataStrctCell = VarDataStrct.fromImgStrctCell(brainImgCell);

[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.kidneys, SUREExperimentCls.windowingKidneys);
kidneyImgCellNov = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortKydneys, mode3D);
kidneyImgCell = ImgStrct.removeNovices(kidneyImgCellNov);
kidneyDataStrctCell = VarDataStrct.fromImgStrctCell(kidneyImgCell);



%%
%Groups of annotators
mainTitle = 'MIN/Mean(Slice Vol) and MAX/Mean(Slice Vol) as function of number of annotators';
annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, 'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
subplot(2,2,1);
VarDataStrct.plotMinMaxVarDiv(liverDataStrctCell,'liver tumors');
subplot(2,2,2);
VarDataStrct.plotMinMaxVarDiv(lungDataStrctCell,'lung tumors');
subplot(2,2,3);
VarDataStrct.plotMinMaxVarDiv(kidneyDataStrctCell,'kidney contour');
subplot(2,2,4);
VarDataStrct.plotMinMaxVarDiv(brainDataStrctCell,'brain hematomas');
    


%%
%Groups of annotators - normalized variability
divByStar = true;
names = {'liver tumors','lung tumors','kindey contour','brain hematomas'};
[varVolCell,varVolMinCell,varVolMaxCell] = VarDataStrct.plotVarDivVarStar({liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},names,divByStar);
for t=1:length(names)
    fprintf('%s ',names{t});
    for z=1:length(varVolCell{t})
        a = num2str(num2str(varVolCell{t}(z),'%.2f'));
        b = num2str(num2str(varVolMaxCell{t}(z)-varVolCell{t}(z),'%.2f'));
        c = num2str(num2str(varVolCell{t}(z)-varVolMinCell{t}(z),'%.2f'));
        fprintf('%s+%s-%s ',a,b,c);
    end
    fprintf('\n',names{t});
end

%%
%Surface distance difference - normalized surface distance
useSurfaceDistance = true;
[varVolCell,varVolMinCell,varVolMaxCell] = VarDataStrct.plotVarDivVarStar({liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},names,divByStar,useSurfaceDistance);
for t=1:length(names)
    fprintf('%s ',names{t});
    for z=1:length(varVolCell{t})
        a = num2str(num2str(varVolCell{t}(z),'%.2f'));
        b = num2str(num2str(varVolMaxCell{t}(z)-varVolCell{t}(z),'%.2f'));
        c = num2str(num2str(varVolCell{t}(z)-varVolMinCell{t}(z),'%.2f'));
        fprintf('%s+%s-%s ',a,b,c);
    end
    fprintf('\n',names{t});
end
%%
%groups of annotators - in one plot
%version 1
params.minMax = false; params.var = true;
mainTitle = 's-vol(variability)/mean annotation s-vol as function of number of annotators';
annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, 'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
subplot(2,2,1);
VarDataStrct.plotMinMaxVarDiv(liverDataStrctCell,'liver tumors',params);
subplot(2,2,2);
VarDataStrct.plotMinMaxVarDiv(lungDataStrctCell,'lung tumors',params);
subplot(2,2,3);
VarDataStrct.plotMinMaxVarDiv(kidneyDataStrctCell,'kidney contour',params);
subplot(2,2,4);
VarDataStrct.plotMinMaxVarDiv(brainDataStrctCell,'brain hematomas',params);

%version 2
VarDataStrct.plotFieldDivValSeveral( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'},...
    'variability');

%%
%%
%Case type and case difficulty
%livers
divByStar = false;
VarDataStrct.plotVarDivVarStar( {liverDataStrctCell(1),liverDataStrctCell(2),liverDataStrctCell(3),liverDataStrctCell(4),liverDataStrctCell(5)}, ...
    {'liver tumor - easy difficulty','liver tumor - easy difficulty','liver tumor - hard difficulty','liver tumor - easy difficulty','liver tumor - medium difficulty'}, ...
    divByStar);
%lungs
VarDataStrct.plotVarDivVarStar({lungDataStrctCell(1),lungDataStrctCell(2),lungDataStrctCell(3),lungDataStrctCell(4),lungDataStrctCell(5)}, ...
    {'lung tumor - hard difficulty','lung tumor - hard difficulty','lung tumor - hard difficulty','lung tumor - medium difficulty','lung tumor - medium difficulty'}, divByStar);


%%
%Surface distance difference - all combined
mainTitle = 'surface distance of possible and consensus, divided by average RECIST';
title(mainTitle);
VarDataStrct.plotFieldDivValSeveral( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'},...
    'surfaceDistanceConPos','meanRECIST');

%%
%Surface distance difference - all in one plot
addText = true; rndStr = '%.2f'; 
annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, 'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
subplot(2,2,1);
title('liver tumors'); ylim([0,0.35]);
[~,minVal,maxVal, meanVal] = VarDataStrct.plotFieldDivVal(liverDataStrctCell,'surfaceDistanceConPos','meanRECIST',addText,rndStr);
 fprintf('liver-tumors ',names{1});VarDataStrct.printMinMaxMeanGaps(meanVal, maxVal, minVal,rndStr);
subplot(2,2,2);
title('lung tumors'); ylim([0,0.35]);
[~,minVal,maxVal, meanVal] = VarDataStrct.plotFieldDivVal(lungDataStrctCell,'surfaceDistanceConPos','meanRECIST',addText,rndStr);
fprintf('lung-tumors ',names{2});VarDataStrct.printMinMaxMeanGaps(meanVal, maxVal, minVal,rndStr);
subplot(2,2,3);
title('kidney contour'); ylim([0,0.35]);
[~,minVal,maxVal, meanVal] = VarDataStrct.plotFieldDivVal(kidneyDataStrctCell,'surfaceDistanceConPos','meanRECIST',addText,rndStr);
fprintf('kidney-contour ',names{3});VarDataStrct.printMinMaxMeanGaps(meanVal, maxVal, minVal,rndStr);
subplot(2,2,4);
title('brain hematomas'); ylim([0,0.35]);
[~,minVal,maxVal, meanVal] = VarDataStrct.plotFieldDivVal(brainDataStrctCell,'surfaceDistanceConPos','meanRECIST',addText,rndStr);
fprintf('brain-hematomas ',names{4});VarDataStrct.printMinMaxMeanGaps(meanVal, maxVal, minVal,rndStr);


%%
%Pairs of annotators - normalized pair overlap
VarDataStrct.plotPairwiseComparison( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    'variability','meanSegSize', ...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'}, ...
    {liverImgCell,lungsImgCell,kidneyImgCell,brainImgCell});
ylim([0,0.7])

%normalized slice volume difference
VarDataStrct.plotPairwiseComparison( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    'volDiff','meanSegSize', ...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'}, ...
    {liverImgCell,lungsImgCell,kidneyImgCell,brainImgCell});
ylim([0,0.7])

VarDataStrct.plotPairwiseComparison( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    'surfaceDistanceConPos','meanRECIST', ...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'}, ...
{liverImgCell,lungsImgCell,kidneyImgCell,brainImgCell});

%Normalized Surface Distance Difference variability
VarDataStrct.plotPairwiseComparison( ...
    {liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},...
    'surfaceDistanceConPos','meanRECIST', ...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'}, ...
{liverImgCell,lungsImgCell,kidneyImgCell,brainImgCell});
ylim([0,0.5])


%%
%Disagreement between annotators - each structure on seperate
%removes annotators from livers and lungs so well have exactly the same
%amount
tempLiverCell = ImgStrct.removeRadiologists(liverImgCell,{'nl1'});%{'nl1','ok'});
tempBrainCell = ImgStrct.removeRadiologists(brainImgCell,{'ns2'});
avgMapLiver = VariabilityExperiment.plotVarDistibution(tempLiverCell,false);
avgMapBrain = VariabilityExperiment.plotVarDistibution(tempBrainCell,false);
avgMapLung = VariabilityExperiment.plotVarDistibution(lungsImgCell,false);
avgMapKidney = VariabilityExperiment.plotVarDistibution(kidneyImgCell,false);
clear tempLiverCell;
clear tempKidneyCell;

% summary of previous run (the average of each structure)
close all;
subplot(2,2,1); VariabilityExperiment.displayMapAsBar(avgMapLiver); title('liver tumors'); xlabel('discrepency level');
subplot(2,2,2); VariabilityExperiment.displayMapAsBar(avgMapLung);title('lung tumors');xlabel('discrepency level');
subplot(2,2,3); VariabilityExperiment.displayMapAsBar(avgMapKidney); title('kidney contours');xlabel('discrepency level');
subplot(2,2,4); VariabilityExperiment.displayMapAsBar(avgMapKidney); title('brain hematomas');xlabel('discrepency level');
mainTitle = 'distribution of number of annotators dissagree with majority within variability (on average)';
annotation('textbox', [0 0.9 1 0.1], 'FontSize',12, ...
    'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');

%%
%%
%compare 2D and 3D data
VarDataStrct.printMinMaxVarForStrcts({liverDataStrctCell,lungDataStrctCell,kidneyDataStrctCell,brainDataStrctCell},N);
VarDataStrct.printMinMaxVarForStrcts({liverDataStrctCell3D,lungDataStrctCell3D,kidneyDataStrctCell3D,brainDataStrctCell3D},N);
N = [8,6,7,7];
mode3D = true;
[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.livers, SUREExperimentCls.windowingLivers);
liverImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongLivers, mode3D);
liverDataStrctCell3D = VarDataStrct.fromImgStrctCell(liverImgCell3D);

[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.lungs, SUREExperimentCls.windowingLungs);
lungImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongLungs, mode3D);
lungDataStrctCell3D = VarDataStrct.fromImgStrctCell(lungImgCell3D);

[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.kidneys, SUREExperimentCls.windowingKidneys);
kidneyImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongKydneys, mode3D);
kidneyDataStrctCell3D = VarDataStrct.fromImgStrctCell(kidneyImgCell3D);

[ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.brains, SUREExperimentCls.windowingBrains);
brainImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongBrains, mode3D);
brainDataStrctCell3D = VarDataStrct.fromImgStrctCell(brainImgCell3D);
%%
%display all annotations of a certain example
im= displaySeperateAnnotations(liverImgCell{3},3,5);


%%
%Case type and case difficulty - each stracture
sortBySeniority = true;

[avgMapLiver, mapLiverCell] = VariabilityExperiment.rankRadiologistsVar(liverImgCellNov,'diceFromOthers',sortBySeniority);
[avgMapLung,mapLungCell] = VariabilityExperiment.rankRadiologistsVar(lungsImgCellNov,'diceFromOthers',sortBySeniority);
[avgMapKidney, mapKidneyCell] = VariabilityExperiment.rankRadiologistsVar(kidneyImgCellNov,'diceFromOthers',sortBySeniority);
[avgMapBrain, mapBrainCell] = VariabilityExperiment.rankRadiologistsVar(brainImgCellNov,'diceFromOthers',sortBySeniority);

VariabilityExperiment.displayMapAsBar(avgMapLiver,VariabilityExperiment.getSeniorityMap(keys(avgMapLiver)),true,true,{'novices','interns','mid-career','experts'});

close all;
mainTitle = 'dice result with mean annotation per Radiologist';
labels = {'novices','interns','mid-career','experts'};
sortBySeniority = true; colorBySeniority = true;
subplot(2,2,1); VariabilityExperiment.displayMapAsBar(avgMapLiver,VariabilityExperiment.getSeniorityMap(keys(avgMapLiver)),sortBySeniority,colorBySeniority,labels); title('liver tumors');
subplot(2,2,2); VariabilityExperiment.displayMapAsBar(avgMapLung,VariabilityExperiment.getSeniorityMap(keys(avgMapLung)),sortBySeniority,colorBySeniority,labels);title('lung tumors');
subplot(2,2,3); VariabilityExperiment.displayMapAsBar(avgMapKidney,VariabilityExperiment.getSeniorityMap(keys(avgMapKidney)),sortBySeniority,colorBySeniority,labels); title('kidney contours');
subplot(2,2,4); VariabilityExperiment.displayMapAsBar(avgMapBrain,VariabilityExperiment.getSeniorityMap(keys(avgMapBrain)),sortBySeniority,colorBySeniority,labels); title('brain hematomas');

annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, ...
    'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');

VariabilityExperiment.calcAverageNumberOfMapShifts(mapLiverCell)
VariabilityExperiment.calcAverageNumberOfMapShifts(mapLungCell)
VariabilityExperiment.calcAverageNumberOfMapShifts(mapKidneyCell)
VariabilityExperiment.calcAverageNumberOfMapShifts(mapBrainCell)
 
tempLiverCell = ImgStrct.removeRadiologists(liverImgCell,{'nl1'});%{'nl1','ok'});
tempLiverCell = ImgStrct.removeRadiologists(tempLiverCell,{'ok'});%{'nl1','ok'});
tempBrainCell = ImgStrct.removeRadiologists(brainImgCell,{'sb2','ns2'}); %%11-9
tempLungCell = ImgStrct.removeRadiologists(lungsImgCell,{'ok'}); %10-9
tempKidnetCell = ImgStrct.removeRadiologists(kidneyImgCell,{'gt'}); %10-9
[~, mapLiverCell] = VariabilityExperiment.rankRadiologistsVar(tempLiverCell,'diceFromOthers',sortBySeniority);
resLiver = VariabilityExperiment.calcMapMean(mapLiverCell);
[~, mapLungCell] = VariabilityExperiment.rankRadiologistsVar(tempLungCell,'diceFromOthers',sortBySeniority);
resLung = VariabilityExperiment.calcMapMean(mapLungCell);
[~, mapKidneyCell] = VariabilityExperiment.rankRadiologistsVar(tempKidnetCell,'diceFromOthers',sortBySeniority);
resKidney = VariabilityExperiment.calcMapMean(mapKidneyCell);
[~, mapBrainCell] = VariabilityExperiment.rankRadiologistsVar(tempBrainCell,'diceFromOthers',sortBySeniority);
resBrain = VariabilityExperiment.calcMapMean(mapBrainCell);


%%
%Average Pearson Correlation
VariabilityExperiment.calcAverageMapPearsonCoeffieicnt(mapLungCell)
VariabilityExperiment.calcAverageMapPearsonCoeffieicnt(mapBrainCell)
VariabilityExperiment.calcAverageMapPearsonCoeffieicnt(mapLiverCell)
VariabilityExperiment.calcAverageMapPearsonCoeffieicnt(mapKidneyCell)

%different measure, instead of pearson correlation
%VariabilityExperiment.calcAverageNumberOfMapShifts(mapLungCell)
%VariabilityExperiment.calcAverageNumberOfMapShifts(mapBrainCell)
%VariabilityExperiment.calcAverageNumberOfMapShifts(mapLiverCell)
%VariabilityExperiment.calcAverageNumberOfMapShifts(mapKidneyCell)

%%
%Case type and case difficulty - per case
close all;
mapAvg = VariabilityExperiment.rankCasesVar({liverImgCell,lungsImgCell,kidneyImgCell,brainImgCell},{'livers','lungs','kidney','brain lesions'});
grMap = containers.Map();
mapKeys = keys(mapAvg);
for z=1:length(mapKeys)
    if strfind(mapKeys{z},'br')
        val = 4;
    elseif strfind(mapKeys{z},'ki')
        val = 3;
    elseif strfind(mapKeys{z},'li')
        val = 1;
    elseif strfind(mapKeys{z},'lu')
        val = 2;
    else
        error('');
    end
    grMap(mapKeys{z}) = val;
end
VariabilityExperiment.displayMapAsBar(mapAvg, grMap, true, true,{'livers','lungs','kidney','brain lesions'});
title('average dice result with mean annotation per case');
xlabel('caes number');
ylabel('average dice result');
%VariabilityExperiment.displayMapAsBar(mapAvg);

%%
sortBySeniority = true;
avgMapLiver = VariabilityExperiment.rankRadiologistsVar(liverImgCellNov,'diceFromOthers',sortBySeniority);
avgMapLung = VariabilityExperiment.rankRadiologistsVar(lungsImgCellNov,'diceFromOthers',sortBySeniority);
%avgMapLung = VariabilityExperiment.rankRadiologistsVar(lungsImgCell([1,2,4,5]),'diceFromOthers',sortBySeniority);
avgMapKidney = VariabilityExperiment.rankRadiologistsVar(kidneyImgCellNov,'diceFromOthers',sortBySeniority);

close all;
mainTitle = 'dice result with mean annotation per Radiologist';
subplot(2,3,1); VariabilityExperiment.displayMapAsBar(avgMapLiver,sortBySeniority); title('liver tumors');
subplot(2,3,2); VariabilityExperiment.displayMapAsBar(avgMapLung,sortBySeniority);title('lung tumors');
subplot(2,3,3); VariabilityExperiment.displayMapAsBar(avgMapKidney,sortBySeniority); title('kidney contours');

annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, ...
    'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
%%
sortBySeniority = false;
avgMapLiver = VariabilityExperiment.rankRadiologistsVar(liverImgCellNov,'diceFromOthers',sortBySeniority);
avgMapLung = VariabilityExperiment.rankRadiologistsVar(lungsImgCell,'diceFromOthers',sortBySeniority);
avgMapKidney = VariabilityExperiment.rankRadiologistsVar(kidneyImgCell,'diceFromOthers',sortBySeniority);

close all;
mainTitle = 'dice result with mean annotation per Radiologist';
subplot(2,3,1); VariabilityExperiment.displayMapAsBar(avgMapLiver); title('liver tumors');
subplot(2,3,2); VariabilityExperiment.displayMapAsBar(avgMapLung);title('lung tumors');
subplot(2,3,3); VariabilityExperiment.displayMapAsBar(avgMapKidney); title('kidney contours');

annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, ...
    'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');

%%

temp = {liverImgCell{4}};
temp{1}.masks = temp{1}.masks(1:6);
tempStrct = VarDataStrct.fromImgStrctCell(temp);
VarDataStrct.plotPairwiseComparison( ...
    {brainImgCell},...
    'surfaceDistanceConPos','meanRECIST', ...
    {'liver tumors','lung tumors','kindey contour','brain hematomas'},...
    {brainImgCell});
VarDataStrct.plotFieldDivValSeveral({liverDataStrctCell},{'liver tumors'}, 'surfaceDistanceConPos','meanRECIST');
res = VarDataStrct.calcPairwiseVarImgCell(liverImgCell{3}, 'surfaceDistanceConPos',true)

Utils.plotLargestSD(brainImgCell{2});



%%
mainTitle = 'var/Mean(Slice Vol)  as function of number of annotators';
annotation('textbox', [0 0.9 1 0.1], 'FontSize',14, 'String', mainTitle, 'EdgeColor', 'none','HorizontalAlignment', 'center');
divByStar = false;

subplot(2,2,1);
[varVolCell,varVolMinCell,varVolMaxCell]  = VarDataStrct.plotVarDivVarStar( {liverDataStrctCell},{'liver tumors'}, divByStar);
subplot(2,2,2);
VarDataStrct.plotVarDivVarStar( {lungDataStrctCell},{'lung tumors'}, divByStar);
subplot(2,2,3);
VarDataStrct.plotVarDivVarStar( {kidneyDataStrctCell},{'kindey contour'}, divByStar);
subplot(2,2,4);
VarDataStrct.plotVarDivVarStar( {brainDataStrctCell},{'brain hematomas'}, divByStar);


