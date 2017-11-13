%filesStrctCell = importdata('matStructs\kidneyLeftSTCroppedNii.mat');
%im = IO.loadFile(filesStrctCell{1}.file);
%seg = IO.loadFile(filesStrctCell{1}.seg);
load('VariabilityEnv.mat');
params = AlgoPriorStrct.getParams();
params.curvature=false; params.intensity=true;
params.shape=true; params.texture=false;
params.shapePriorAtlas = ShapeAtlas.LeftCloseDefaultAtlas();
params.shapePriorMask = prior.shapePriorMask;
params.min=true; params.max=false; params.mean = false;
%algoPrior = Prior.calculatePriors(im,seg,params);



l = 60;
I = im.img(:,:,l);
S = logical(seg.img(:,:,l));
sliceParams = params;
sliceParams.shapePriorMask=params.shapePriorMask(:,:,l);
func = @(segmentationIn, roi) VariabilityEstimator.calculatePriorsWrapper(I,segmentationIn,sliceParams,'min',roi);
priorDisp = Utils.displayCertaintyUncertainty2_3D(I,prior.min(:,:,l),Utils.getBoundries(S,false));

estimator = VariabilityEstimator(I,S,func,0.5);
[uncertaintyRegions,mask] = estimator.extractUncertaintyRegions(prior.min(:,:,l));
[res] = estimator.measureRegionSensitivity(uncertaintyRegions.PixelIdxList{2},uncertaintyRegions.angle{2});
[variabilityData] = estimator.measureRegionsSensitivity(uncertaintyRegions);
variabilityData = estimator.evaluateVariability(variabilityData);

varMask2D = estimator.getVariabilityMask(variabilityData);
varMask2D = imclose(varMask2D,strel('disk',3));
res = Utils.displaySegmentation(priorDisp,varMask2D,[0,0,1]);
res = Utils.displaySegmentation(estimator.im,estimator.seg,[1,0,0]);
res = Utils.displayCertaintyUncertainty2_3D(res,varMask2D,Utils.getBoundries(varMask2D,true));

segRes = Utils.displaySegmentation(I,S,[0,0,1]);
imshow(squeeze(segRes));
segMask

IO.imshow3D(res)

res = estimator.getSensitivityDbgRes(variabilityData); IO.imshow3D(res);
estimator.displayUncertaintyRegions(squeeze(priorDisp),variabilityData);
res = estimator.estimate(im.img(:,:,60), seg.img(:,:,60));
[lines, dirs, simplfyContour] = estimator.getContourPatches(seg.img(:,:,60));

drorSeg = IO.loadFile('data\simpleHypoCheck\drorSegs_10000101_1_29663_305_Dror.nii.gz');
uncertaintyMask = xor(drorSeg.img>0,seg.img>0);

params.shapePriorMask = prior.shapePriorMask;
[varMask, varData, dbgSensRes] = VariabilityEstimator.evaluate3DVarMask(im.img,logical(seg.img),prior.min,params,0.5);
varMask = imclose(varMask,strel('disk',3));
[priorDisp,overlay] = Utils.displayCertaintyUncertainty2_3D(im.img,prior.min,Utils.getBoundries(seg.img,false));
res = Utils.displaySegmentation(priorDisp,varMask,[0,3/5,1],false);
res = Utils.displayCertaintyUncertainty2_3D(res,varMask,Utils.getBoundries(varMask,true));


figure,IO.imshow3D(res)

zFirst = find(sum(sum(seg.img,1),2) > 0,1,'first');
zLast = find(sum(sum(seg.img,1),2) > 0,1,'last');
IO.saveAsTiff('win20.tiff',res(:,:,zFirst:zLast,:));
figure,IO.imshow3D(dbgSensRes)

res = Utils.displaySegmentation(priorDisp,uncertaintyMask,[0,3/5,1],true);
IO.saveAsTiff('actualVar.tif',res(:,:,zFirst:zLast,:));

out = Utils.getContourLineOverlay(prior.min(:,:,60),im.img(:,:,60));

[priorDisp,overlay] = Utils.displayCertaintyUncertainty2_3D(im.img,prior.min,Utils.getBoundries(seg.img,false));
lineStrctCell = LineStruct.getLineStructCell(overlay);
lineStrctCell = LineStruct.getLineStructCell(varMask,[0,0,1],lineStrctCell);
LineDisplay.ui(im.img,lineStrctCell);

lineStrctCell2 = LineStruct.getLineStructCell(overlay);
lineStrctCell2 = LineStruct.getLineStructCell(uncertaintyMask,[0,0,1],lineStrctCell2);
LineDisplay.ui(im.img,lineStrctCell2);

lineStrctCell3 = LineStruct.getLineStructCell(drorSeg.img,[0,0,1]);
lineStrctCell3 = LineStruct.getLineStructCell(seg.img,[0,1,1],lineStrctCell3);
LineDisplay.ui(im.img,lineStrctCell3);

out = VariabilityEstimator.calcMeanShape(drorSeg.img,seg.img);