basefolder = '/cs/casmip/clara.herscu/git/thesis/SUREExperimentDataForUpload/';
annotationDirs = strcat([basefolder, 'anottations/'], SUREExperimentCls.novicesDirs);
dataDir = [basefolder, 'data/'];
mode3D = false;

[ dbFiles ] = generateDbStrctFromFolder(annotationDirs, dataDir, SUREExperimentCls.lungs, SUREExperimentCls.windowingLungs);
imgCell = ImgStrct.dbFileStrctToImgStrctCell(dbFiles, SUREExperimentCls.slicesShortLungs, mode3D);

imNum = 1;
im = imgCell{imNum}.img;
segs = imgCell{imNum}.masks;
meanSeg = stapleWrapper(segs) > 0.5;

% extracting only boundary points
segBoundary = Utils.getBoundries(meanSeg);
% getting their coordinates
[row, col] = find(segBoundary);

% this doesn't work
segPoints = [col(:) row(:)];

% % displaying the points for debug purposes
% imshow(im)
% hold on
% plot(segPoints(:,2), segPoints(:,1), 'b.')

Options = struct;
Options.Iterations = 10;
Options.nPoints = 70;
Options.Delta = 1;
Option.Verbose = true;
[O,J]=Snake2D(im,segPoints_subsampled,Options);

% view results
LineDisplay.displayMasks(im, J);
