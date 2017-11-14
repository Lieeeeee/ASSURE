basefolder = '/cs/casmip/clara.herscu/git/thesis/SUREExperimentDataForUpload/';
annotationDirs = strcat([basefolder, 'anottations/'], SUREExperimentCls.novicesDirs);
dataDir = [basefolder, 'data/'];
%set to true to move to 3D mode
mode3D = false;
%reads lungs tumore files
[ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.lungs, SUREExperimentCls.windowingLungs);
imgCell = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortLungs, mode3D);

% [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.livers, SUREExperimentCls.windowingLivers);
% imgCell = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongLivers, mode3D);
% 
% [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.kidneys, SUREExperimentCls.windowingKidneys);
% imgCell = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongKydneys, mode3D);
% 
% [ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.brains, SUREExperimentCls.windowingBrains);
% imgCell3D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongBrains, mode3D);
% imgCell2D = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.sliceLongBrains, false);

% [ dbFiles ] = generateDbStrctFromFolder(annotationsDirs, dataDir, SUREExperimentCls.kidneys, SUREExperimentCls.windowingKidneys);
% imgCell = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortKydneys, mode3D);

% for imNum = 1:size(imgCell, 2)
imNum = 2;
%     %shows actual variability
%     Utils.showVariability(imgCell{imNum})
%     %shows annotation
%     Utils.showAnnotations(imgCell{imNum})
%     %shows possible and consensus
%     Utils.showConsensusPossibleAsLines(imgCell{imNum})

%%
%initalizes params
params = AlgoPriorStrct.getParams();
params.curvature=true; 
%intensity filter parameters
params.intensity=false; %intensity prior where the thrshold is global
params.intensityLocal = true; %intensity prior where the threshold is local
params.filterSize = 1;
%parms.intensityLocal = ... %optional - set a threshold for intensity local prior 

params.texture=false;

%shape atlas parameters
params.shape=false; 
%params.shapePriorAtlas = ShapeAtlas.LeftCloseDefaultAtlas();
%integrator function
params.min = true;

%%
%perform variability and priors calculation

im = imgCell{imNum}.img;
% perform histogram equalization
% [im, ~] = histeq(im);
% seg = imgCell{imNum}.masks{1};
% varMask = VariabilityEstimator.evaluate3DVarMaskTrivial(im,seg,prior.min,params,0.5);
segs = imgCell{imNum}.masks;
% meanSeg = VariabilityEstimator.calcMeanShapeMultSegs(segs);
meanSeg = stapleWrapper(segs) > 0.5;

%calculate all priors according to params
prior = Prior.calculatePriors(im,meanSeg,params);
%calculate just intensity prior
intensityPrior = Prior.intensityPriorLocal(im,seg,2);
%     %display priors result
%     LineDisplay.displayQuality(im,meanSeg,prior.curvature)
%     LineDisplay.displayQuality(im,meanSeg,prior.intensityLocal)
%     LineDisplay.displayQuality(im,meanSeg,prior.min)

% calculate variability from mean segmentation
varMask_2 = VariabilityEstimator.evaluate3DVarMask(im,meanSeg,prior.min,params,0.5); % dror's way
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(im,meanSeg); % my way

%display the variability result
% LineDisplay.displayVariability(im,meanSeg,prior.min,varMask);
LineDisplay.displayVariabilityWithoutSeg(im, varMask_2); % dror's way
LineDisplay.displayVariabilityWithoutSeg(im, varMask); % my way

% display the true variability
[~,~,trueVarMask] = Utils.calcUnionIntersection(imgCell{imNum}.masks);
LineDisplay.displayVariabilityWithoutSeg(im, trueVarMask);

%% measuring success
% normalized volume overlap = 2|est & gt| / |est| + |gt|

% percentage diameter error = (mean_diameter(est)-mean_diameter(gt))/mean_diameter(gt)

% percentage volume error = (vol(est) - vol(gt))/vol(gt)

% end