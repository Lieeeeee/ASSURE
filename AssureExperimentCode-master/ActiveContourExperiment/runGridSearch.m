% run a grid search on every model we want to find the parameters of
LOAD_DATA = true;
DISPLAY = false;
WINDOWS = false;
DEBUG = false;

% load data
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

experiment_folder = '/cs/casmip/clara.herscu/git/thesis/figs/snakes_jan18/experiment_010118/';

% lung
kernelSize = 2; Tlength = 3;
girdSearchSnake( lungImgCell3D, [experiment_folder 'lung/'], kernelSize, Tlength );

% liver
kernelSize = 4; Tlength = 3;
girdSearchSnake( liverImgCell3D, [experiment_folder 'liver/'], kernelSize, Tlength );

% kidney 
kernelSize = 3; Tlength = 7;
girdSearchSnake( kidneyImgCell3D, [experiment_folder 'kidney/'], kernelSize, Tlength );

% brain
kernelSize = 3; Tlength = 5;
girdSearchSnake( brainImgCell3D, [experiment_folder 'brain/'], kernelSize, Tlength );
