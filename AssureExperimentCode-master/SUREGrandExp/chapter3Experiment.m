LOAD_DATA = true;
DISPLAY = false;


if LOAD_DATA
    load('chapter3Env.mat');
else 
    annotationDirs = {'C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\largeAmount'}; %...
    %'C:\Users\drorcohe\Desktop\SUREGrandExperiment\anottations\smallAmount';
    dataDir = 'C:\Users\drorcohe\Desktop\SUREGrandExperiment\data';
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

%% initializing with mean segmentation
% kernelSize = 2; Tlength = 3;
% [resLung, outMasksLung] =  VariabilityExperiment.holdExperiment(lungImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/dror_lung/',Tlength);%,'lung'); new
% printRes(resLung)
% 
% kernelSize = 4; Tlength = 3;
% resLiver = VariabilityExperiment.holdExperiment(liverImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/dror_liver/',Tlength);
% printRes(resLiver)
% 
% kernelSize = 3; Tlength = 7;
% [resKidney, outMasksKidney] =  VariabilityExperiment.holdExperiment(kidneyImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/dror_kidney/',Tlength);
% 
% printRes(resKidney)
% 
% kernelSize = 3; Tlength = 5;
% [resBrain, outMasksBrain] = VariabilityExperiment.holdExperiment(brainImgCell3D,kernelSize,[],'/cs/casmip/clara.herscu/git/thesis/figs/activeContour_nov17/experiment_131117/dror_brain/',Tlength); %new
% printRes(resBrain)
% 
% 
% % save('estimatedVariabilityEnv.mat')

%% initializing with different segmentation every time
experiment_folder = '/cs/casmip/clara.herscu/git/thesis/figs/activeContour_dec17/experiment_061217/dror_different_init/';
allResStruct = struct;

kernelSize = 1; Tlength = 3;
for segNum = [1,3:7]
    [resLung, outMasksLung] =  VariabilityExperiment.holdExperiment(lungImgCell3D,kernelSize,[],[experiment_folder 'lung/' num2str(segNum) '/'],Tlength,0,0,segNum);
    printRes(resLung)
    allResStruct.lungRes{segNum} = resLung;
    allResStruct.lungOutMasks{segNum} = outMasksLung;
end


kernelSize = 4; Tlength = 3;
for segNum = 1:8
    [resLiver, outMasksLiver] = VariabilityExperiment.holdExperiment(liverImgCell3D,kernelSize,[],[experiment_folder 'liver/' num2str(segNum) '/'],Tlength,0,0,segNum);
    printRes(resLiver)
    allResStruct.liverRes{segNum} = resLiver;
    allResStruct.liverOutMasks{segNum} = outMasksLiver;
end


kernelSize = 3; Tlength = 7;
for segNum = 1:7
    [resKidney, outMasksKidney] =  VariabilityExperiment.holdExperiment(kidneyImgCell3D,kernelSize,[],[experiment_folder 'kidney/' num2str(segNum) '/'],Tlength,0,0,segNum);
    printRes(resKidney)
    allResStruct.kidneyRes{segNum} = resKidney;
    allResStruct.kidneyOutMasks{segNum} = outMasksKidney;
end

kernelSize = 3; Tlength = 5;
for segNum = 1:7
    [resBrain, outMasksBrain] = VariabilityExperiment.holdExperiment(brainImgCell3D,kernelSize,[],[experiment_folder 'brain/' num2str(segNum) '/'],Tlength,0,0,segNum);
    printRes(resBrain)
    allResStruct.brainRes{segNum} = resBrain;
    allResStruct.brainOutMasks{segNum} = outMasksBrain;
end

save('AssureExperimentCode-master/Dror_DifferentInitializations_Experiment061217.mat')
