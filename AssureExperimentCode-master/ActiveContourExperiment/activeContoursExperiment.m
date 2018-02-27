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

%% using mean segmentations
experiment_folder = '/cs/casmip/clara.herscu/git/thesis/figs/activeContour_feb18/experiment_270218/';
% [OptionsIn, OptionsOut] = ActiveContourOptions.getLungOptions();
params = [10, 0, 0, 10, 0, 0]; 
iterIn = params(1); cIn = params(2); sIn = params(3); 
OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);
iterOut = params(4); cOut = params(5); sOut = params(6); 
OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);

kernelSize = 2; Tlength = 3;
[resLung, outMasksLung] =  VariabilityExperiment.holdExperiment(lungImgCell3D,kernelSize,[],[experiment_folder 'lung/'],Tlength,0,1,0,OptionsIn, OptionsOut);
printRes(resLung)
% 
% kernelSize = 4; Tlength = 3;
% [OptionsIn, OptionsOut] = ActiveContourOptions.getLiverOptions();
% resLiver = VariabilityExperiment.holdExperiment(liverImgCell3D,kernelSize,[],[experiment_folder 'liver/'],Tlength,0,1,0,OptionsIn, OptionsOut);
% printRes(resLiver)
% 
% kernelSize = 3; Tlength = 7;
% [OptionsIn, OptionsOut] = ActiveContourOptions.getKidneyOptions();
% [resKidney, outMasksKidney] =  VariabilityExperiment.holdExperiment(kidneyImgCell3D,kernelSize,[],[experiment_folder 'kidney/'],Tlength,0,1,0,OptionsIn, OptionsOut);
% printRes(resKidney)
% 
% kernelSize = 3; Tlength = 5;
% [OptionsIn, OptionsOut] = ActiveContourOptions.getBrainOptions();
% [resBrain, outMasksBrain] = VariabilityExperiment.holdExperiment(brainImgCell3D,kernelSize,[],[experiment_folder 'brain/'],Tlength,0,1,0,OptionsIn, OptionsOut);
% printRes(resBrain)
% 
% n_iterations = 4; % 5
% contraction_param_out = -0.5; % -0.5
% smooth_param_out = 0.1; % 0
% 
% save('AssureExperimentCode-master/ActiveContourExperiment/activeContourEstimationExperiment_171217_meanSegInit.mat')

%% using other segmentation every time
% experiment_folder = '/cs/casmip/clara.herscu/git/thesis/figs/activeContour_dec17/experiment_171217_with_smoothing/different_init/';
% allResStruct = struct;
% 
% kernelSize = 1; Tlength = 3;
% for segNum = [1, 3:7]
%     [resLung, outMasksLung] =  VariabilityExperiment.holdExperiment(lungImgCell3D,kernelSize,[],[experiment_folder 'lung/' num2str(segNum) '/'],Tlength,0,1,segNum);
%     printRes(resLung)
%     allResStruct.lungRes{segNum} = resLung;
%     allResStruct.lungOutMasks{segNum} = outMasksLung;
% end
% 
% 
% kernelSize = 4; Tlength = 3;
% for segNum = 1:8
%     [resLiver, outMasksLiver] = VariabilityExperiment.holdExperiment(liverImgCell3D,kernelSize,[],[experiment_folder 'liver/' num2str(segNum) '/'],Tlength,0,1,segNum);
%     printRes(resLiver)
%     allResStruct.liverRes{segNum} = resLiver;
%     allResStruct.liverOutMasks{segNum} = outMasksLiver;
% end
% 
% 
% kernelSize = 3; Tlength = 7;
% for segNum = 1:7
%     [resKidney, outMasksKidney] =  VariabilityExperiment.holdExperiment(kidneyImgCell3D,kernelSize,[],[experiment_folder 'kidney/' num2str(segNum) '/'],Tlength,0,1,segNum);
%     printRes(resKidney)
%     allResStruct.kidneyRes{segNum} = resKidney;
%     allResStruct.kidneyOutMasks{segNum} = outMasksKidney;
% end
% 
% kernelSize = 3; Tlength = 5;
% for segNum = 1:7
%     [resBrain, outMasksBrain] = VariabilityExperiment.holdExperiment(brainImgCell3D,kernelSize,[],[experiment_folder 'brain/' num2str(segNum) '/'],Tlength,0,1,segNum);
%     printRes(resBrain)
%     allResStruct.brainRes{segNum} = resBrain;
%     allResStruct.brainOutMasks{segNum} = outMasksBrain;
% end
% 
% % also saving active contour parameter
% n_iterations = 4;
% contraction_param_out = -0.5;
% smooth_param_out = 0.1;
% 
% save('AssureExperimentCode-master/ActiveContourExperiment/activeContourDifferentInitializationsExperiment171217_with_smoothing_.mat')