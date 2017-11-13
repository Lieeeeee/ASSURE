%%shows annotations
Utils.displaySeveralStructs(@Utils.showIm, {liverImgCell{4},lungsImgCell{3},kidneyImgCell{6},brainImgCell{1}},[3,4,3,5]);
Utils.displaySeveralStructs(@Utils.showConsensusPossibleAsLines, {liverImgCell{4},lungsImgCell{3},kidneyImgCell{6},brainImgCell{1}},[3,4,3,5]);
Utils.displaySeveralStructs(@Utils.showAnnotations, {liverImgCell{4},lungsImgCell{3},kidneyImgCell{6},brainImgCell{1}},[3,4,3,5]);
Utils.displaySeveralStructs(@Utils.showVariability, {liverImgCell{4},lungsImgCell{3},kidneyImgCell{6},brainImgCell{1}},[3,4,3,5]);

Utils.displaySeveralStructs(@Utils.showVariability, {liverImgCell{5},liverImgCell{4},lungsImgCell{1},lungsImgCell{5}},[4,3,5,5]);
Utils.displaySeveralStructs(@Utils.showIm, {liverImgCell{5},liverImgCell{4},lungsImgCell{1},lungsImgCell{5}},[4,3,5,5]);

Utils.displaySeveralStructs(@Utils.showPossibleConsensusSurfaceDistance, {liverImgCell{4},lungsImgCell{3},kidneyImgCell{6},brainImgCell{2}},[3,4,3,4]);


Utils.displaySeveralStructs(@Utils.showAnnotations,{kidneyImgCell{3},brainImgCell{1}},[10,5]);
LineDisplay.displayMasks(imgStrct.img,{intersection});



[lowestLiverVar, lowestLiverVarVal] = Utils.findLowestVarSlice(liverImgCell, true);
[lowestLiverSD, lowestLiverSDVal] = Utils.findLowestVarSlice(liverImgCell, false);
[lowestLungVar, lowestLungVarVal] = Utils.findLowestVarSlice(lungsImgCell, true);
[lowestLungSD, lowestLungSDVal] = Utils.findLowestVarSlice(lungsImgCell, false);
[lowestKidneyVar, lowestKidneyVarVal] = Utils.findLowestVarSlice(kidneyImgCell, true);
lowestKidneyVar = ImgStrct.cropImgStrct(lowestKidneyVar,1,183,200,363);
[lowestKidneySD, lowestKidneySDVal] = Utils.findLowestVarSlice(kidneyImgCell, false);
[lowestBrainVar, lowestBrainVarVal] = Utils.findLowestVarSlice(brainImgCell, true);
[lowestBrainSD, lowestBrainSDVal] = Utils.findLowestVarSlice(brainImgCell, false);
Utils.displaySeveralStructs(@Utils.showIm, {lowestLiverVar,lowestLungVar,lowestKidneyVar,lowestBrainVar},[1,1,1,1]);
labels = {num2str(lowestLiverVarVal),num2str(lowestLungVarVal),num2str(lowestKidneyVarVal),num2str(lowestBrainVarVal)};
Utils.displaySeveralStructs(@Utils.showVariability, {lowestLiverVar,lowestLungVar,lowestKidneyVar,lowestBrainVar},[1,1,1,1],labels);
labels = {num2str(lowestLiverSDVal),num2str(lowestLungSDVal),num2str(lowestKidneySDVal),num2str(lowestBrainSDVal)};
Utils.displaySeveralStructs(@Utils.showPossibleConsensusSurfaceDistance, {lowestLiverSD,lowestLungSD,lowestKidneySD,lowestBrainSD},[1,1,1,1],labels);

Utils.showVariability(lowestKidneyVar)
Utils.showIm(lowestKidneyVar)
Utils.showAnnotations(lowestKidneyVar)
plotMin =false; Utils.plotLargestSD(lowestKidneyVar, 2, plotMin);
plotMin =true; Utils.plotLargestSD(lowestKidneyVar, 2, plotMin);
Utils.showPossibleConsensusSurfaceDistance(lowestKidneySD)
Utils.showIm(lowestKidneySD)
Utils.showAnnotations(lowestKidneySD)
plotMin =false; Utils.plotLargestSD(lowestKidneySD, 2, plotMin);
plotMin =true; Utils.plotLargestSD(lowestKidneySD, 2, plotMin);
Utils.showPossibleConsensusSurfaceDistance(lowestKidneySD)
im0 = LineDisplay.getCroppedFrameFromFigure();

[possible,consensus,~] = Utils.calcUnionIntersection(lowestKidneySD.masks);
[SD,meanSD] = Utils.surfaceDistance(possible,consensus);

close all;


