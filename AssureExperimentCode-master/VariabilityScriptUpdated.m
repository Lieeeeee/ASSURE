%filesStrctCell = importdata('matstructs\kidneyLeftGT.mat');
%for tt=1:length(filesStrctCell)
%filesStrctCell{tt}.algoFiles=[];
%end
%filesStrctCellToEval = IO.copyFileStrctToDir(filesStrctCell, 'data\varEstimation\');
%rmvList = [];
%for tt=1:length(filesStrctCellToEval)
%    if isempty(filesStrctCellToEval{tt}.seg)
%       rmvList = [rmvList,tt];
%    end
%end
%filesStrctCellToEval(rmvList) = [];
%IO.cropFileStrct(filesStrctCellToEval);

%filesStrctCell = importdata('data\varEstimation\kidneysLeftSTCroppedNii.mat');
filesStrctCell = importdata('data\varEstimation\kidneysLeftGTCroppedNii.mat');

filesStrctCell(3) = [];
params = AlgoPriorStrct.getParams();
params.curvature=false; params.intensity=true;
params.shape=false; params.texture=false;
params.shapePriorAtlas = ShapeAtlas.LeftCloseDefaultAtlas();
%params.shapePriorMask = prior.shapePriorMask;
params.min=true; params.max=false; params.mean = false;
gap = 5;


im1 = IO.loadFile(filesStrctCell{1}.file);
seg1 = IO.loadFile(filesStrctCell{1}.seg);





NSamples = 20;


samplesVar2D = ones(NSamples,150)*-1;
samplesVar = ones(NSamples,1)*-1;

saveFile = false;
for tt=1:NSamples
    tt
    im = IO.loadFile(filesStrctCell{tt}.file);
    seg = IO.loadFile(filesStrctCell{tt}.seg);
    
    continueFlag = false;
    for z=1:size(seg.img,3)
        seg.img(:,:,z) = imfill(seg.img(:,:,z),'holes');
    end
    if continueFlag
        continue;
    end
    
    prior = Prior.calculatePriors(im,seg,params);
    
    
    varMask = VariabilityEstimator.evaluate3DVarMask(im.img,seg.img,prior.min,params,0.5);
    
    zStart = find(sum(sum(seg.img,1),2)>0,1,'first') + gap;
    zEnd = find(sum(sum(seg.img,1),2)>0,1,'last') - gap;
    
    varMaskInRoi = varMask(:,:,zStart:zEnd);
    segInRoi = seg.img(:,:,zStart:zEnd);
    varRatio = sum(sum(varMaskInRoi,1),2); %sum(sum(segInRoi,1),2);
    samplesVar2D(tt,zStart:zEnd) = varRatio;
    varRatio3D = sum(varMaskInRoi(:))/size(segInRoi,3);%sum(segInRoi(:));
    samplesVar(tt) = varRatio3D;
    samplesVar
    if saveFile
        tiffName = ['data\varEstimation\out\' num2str(tt) '.tiff'];
        try
            LineDisplay.saveVariabilityAsTiff(im.img(:,:,zStart:zEnd),seg.img(:,:,zStart:zEnd),prior.min(:,:,zStart:zEnd),varMask(:,:,zStart:zEnd),tiffName);
        catch
            LineDisplay.saveVariabilityAsTiff(im.img(:,:,zStart:zEnd),seg.img(:,:,zStart:zEnd),prior.min(:,:,zStart:zEnd),varMask(:,:,zStart:zEnd),tiffName);
        end
    end
    %   LineDisplay.displayVariability(im.img(:,:,zStart:zEnd),seg.img(:,:,zStart:zEnd),prior.min(:,:,zStart:zEnd),varMask(:,:,zStart:zEnd));
    
    save('data\varEstimation\out\sampleVar.mat','samplesVar');
    save('data\varEstimation\out\samplesVar2D.mat','samplesVar2D');
end


 varMask2D = VariabilityEstimator.evaluate3DVarMask(im.img(:,:,65),seg.img(:,:,65),prior.min(:,:,65),params,0.5);
  LineDisplay.displayVariability(im.img(:,:,65),seg.img(:,:,65),prior.min(:,:,65),varMask2D);
     

%filesStrctCell = importdata('matStructs\kidneyLeftSTCroppedNii.mat');
%im = IO.loadFile(filesStrctCell{1}.file);
%seg = IO.loadFile(filesStrctCell{1}.seg);
%load('VariabilityEnv.mat');
params = AlgoPriorStrct.getParams();
params.curvature=false; params.intensity=true;
params.shape=false; params.texture=false;
params.shapePriorAtlas = ShapeAtlas.LeftCloseDefaultAtlas();
params.shapePriorMask = prior.shapePriorMask;
params.min=true; params.max=false; params.mean = false;
%algoPrior = Prior.calculatePriors(im,seg,params);

drorSeg = IO.loadFile('data\simpleHypoCheck\drorSegs_10000101_1_29663_305_Dror.nii.gz');

%claculates mean seg
meanSeg = VariabilityEstimator.calcMeanShape(drorSeg.img,seg.img);
meanSegPrior = Prior.calculatePriors(im.img,meanSeg,params);

%calc variability from mean seg
[varMask, varData, dbgSensRes] = VariabilityEstimator.evaluate3DVarMask(im.img,meanSeg,meanSegPrior.min,params,0.5);
LineDisplay.displayVariability(im.img,meanSeg,meanSegPrior.min,varMask);

%display variability on orig seg
LineDisplay.displayVariability(im.img,seg.img,prior.min,varMask);

%calc variability on orig seg
varMaskSeg= VariabilityEstimator.evaluate3DVarMask(im.img,seg.img,prior.min,params,0.5);
LineDisplay.displayVariability(im.img,seg.img,prior.min,varMaskSeg);

%shows uncertainty

uncertaintyMask = xor(drorSeg.img>0,seg.img>0);
lineStrctCell = LineStruct.getLineStructCell(drorSeg.img,[0,0,1]);
lineStrctCell = LineStruct.getLineStructCell(seg.img,[0,1,1],lineStrctCell);
LineDisplay.ui(im.img,lineStrctCell);

