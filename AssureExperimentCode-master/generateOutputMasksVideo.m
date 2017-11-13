
load('estimatedVariabilityEnv.mat');
load('chapter3Env.mat');

dataDir = 'C:\Users\drorcohe\Desktop\SUREGrandExperiment\data';
caseType= 4;
saveAsVid = true;
saveEstimated = true;
saveActual = true;
%caseNum = 1;

close all
currentWindowing = SUREExperimentCls.windowingLivers;
currentImgCell = liverImgCell3D;
slicesLong = SUREExperimentCls.sliceLongLivers;
outMask = outMasksLiver;
if caseType==2
    currentImgCell = lungImgCell3D;
    slicesLong = SUREExperimentCls.sliceLongLungs;
    outMask = outMasksLung;
elseif caseType==3
    currentImgCell = kidneyImgCell3D;
    slicesLong = SUREExperimentCls.sliceLongKydneys;
    outMask = outMasksKidney;
elseif caseType==4
    currentImgCell = brainImgCell3D;
    outMask = outMasksBrain;
    slicesLong = SUREExperimentCls.sliceLongBrains;
end

for caseNum=1:length(currentImgCell)
    vol = currentImgCell{caseNum}.img;
    
    
    [~,b,~] = fileparts(currentImgCell{caseNum}.imgFileName);
    [lineStrctCell, colors] = LineDisplay.displaySegsOverlay(vol,currentImgCell{caseNum}.masks,[],false);
    if saveEstimated
        volMask =outMask{caseNum};
        vid1Name = [ b 'variabilityMaskEstimated.avi'];
        LineDisplay.saveLineStrctCellAsTiff(vol, lineStrctCell,vid1Name ,saveAsVid,volMask)
    end
    if saveActual
        vid2Name = [ b 'variabilityMaskActual.avi'];
        [~,~,volMaskActual] = Utils.calcUnionIntersection(currentImgCell{caseNum}.masks);
        LineDisplay.saveLineStrctCellAsTiff(vol, lineStrctCell, vid2Name ,saveAsVid,volMaskActual)
    end
    IO.saveAsVid([ currentFileStrct{caseNum} 'img.avi'],vol);
    if saveActual && saveEstimated
        IO.uniteVideos(vid1Name,vid2Name,[ b 'bothVar.avi']);
    end
end
