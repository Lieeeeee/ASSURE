LOAD_DATA = true;
DISPLAY = false;


if LOAD_DATA
    load('chapter3Env.mat');
else
    basefolder = '/cs/casmip/clara.herscu/git/thesis/SUREExperimentDataForUpload/';
    annotationsDirs = strcat([basefolder, 'anottations/'], SUREExperimentCls.novicesDirs);
    dataDir = [basefolder, 'data/'];
%     annotationDirs = {'C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\largeAmount'}; %...
    %'C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\smallAmount';
%     dataDir = 'C:\Users\drorcohe\Desktop\SUREGrandExperiment\data';
    mode3D = true;
    
    [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.livers, SUREExperimentCls.windowingLivers);
    liverImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongLivers, mode3D);
    
    [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.kidneys, SUREExperimentCls.windowingKidneys);
    kidneyImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongKydneys, mode3D);
    
    [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.lungs, SUREExperimentCls.windowingLungs);
    lungImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongLungs,mode3D);
    
    [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.brains, SUREExperimentCls.windowingBrains);
    brainImgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongBrains, mode3D);
    brainImgCell2D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongBrains, false);
end

kernelSize = 2; Tlength = 3;
[resLung, outMasksLung] =  VariabilityExperiment.holdExperiment(lungImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/lung/',Tlength,0,1);%,'lung'); new
printRes(resLung)

kernelSize = 4; Tlength = 3;
resLiver = VariabilityExperiment.holdExperiment(liverImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/liver/',Tlength,0,1);
printRes(resLiver)

kernelSize = 3; Tlength = 7;
[resKidney, outMasksKidney] =  VariabilityExperiment.holdExperiment(kidneyImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/kidney/',Tlength,0,1);
printRes(resKidney)

kernelSize = 3; Tlength = 5;
[resBrain, outMasksBrain] = VariabilityExperiment.holdExperiment(brainImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/brain/',Tlength,0,1); %new
printRes(resBrain)


save('AssureExperimentCode-master/activeContourEstimationExperiment.mat')




