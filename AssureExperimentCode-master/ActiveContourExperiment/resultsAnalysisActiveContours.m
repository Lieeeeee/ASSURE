%% mean segmentation initialization
% load('/cs/casmip/clara.herscu/git/thesis/AssureExperimentCode-master/ActiveContourExperiment/activeContourEstimationExperiment_141217_meanSegInit.mat')
load('/cs/casmip/clara.herscu/git/thesis/AssureExperimentCode-master/ActiveContourExperiment/activeContourEstimationExperiment_171217_meanSegInit.mat')

display('========== lung ==========');
printRes(resLung);

display('========== liver ==========');
printRes(resLiver);

display('========== kidney ==========');
printRes(resKidney);

display('========== brain ==========');
printRes(resBrain);

clear;
%% different initializations
load('/cs/casmip/clara.herscu/git/thesis/AssureExperimentCode-master/ActiveContourExperiment/activeContourDifferentInitializationsExperiment141217_with_smoothing_.mat')

display('========== lung ==========');
for segNum = [1, 3:7]
    display(['initialized with seg no. ' num2str(segNum)]);
    printRes(allResStruct.lungRes{segNum});
end

display('========== liver ==========');
for segNum = 1:8
    display(['initialized with seg no. ' num2str(segNum)]);
    printRes(allResStruct.liverRes{segNum});
end

display('========== kidney ==========');
for segNum = 1:7
    display(['initialized with seg no. ' num2str(segNum)]);
    printRes(allResStruct.kidneyRes{segNum});
end

display('========== brain ==========');
for segNum = 1:7
    display(['initialized with seg no. ' num2str(segNum)]);
    printRes(allResStruct.brainRes{segNum});
end