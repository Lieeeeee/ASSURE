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

% segPoints = [col(:) row(:)];

% % displaying the points for debug purposes
% imshow(im)
% hold on
% plot(segPoints(:,2), segPoints(:,1), 'b.')

%% % implementation 1 -- works, but contour points have to be ordered clockwise
% order points
% [row_clockwise, col_clockwise] = poly2cw(row_clockwise, col_clockwise);
colCenter = mean(col);
rowCenter = mean(row);
angles = atan2d((row-rowCenter), (col-colCenter));
[~, sortedIndexes] = sort(angles);
col_clockwise = col(sortedIndexes);  % Reorder x and y with the new sort order.
row_clockwise = row(sortedIndexes);



segPoints = [col_clockwise(:) row_clockwise(:)];

% displaying initial contour
imshow(im);
hold on
plot(segPoints(:,1), segPoints(:,2), 'b');
hold on 
plot(rowCenter, colCenter, 'g+', 'MarkerSize', 15);

% setting the options
Options = struct;
Options.Iterations = 1;
Options.nPoints = 94;
Options.Delta = 1;
Options.Verbose = true;
Options.Wline = -0.04;
Options.Kappa = -2;
[O,J]=Snake2D(im,fliplr(segPoints),Options);

% view results
LineDisplay.displayMasks(im, J);
%% % implementation 2 -- doesn't seem to work
% % image: This is the image data
% % xs, ys: The initial snake coordinates
% % alpha: Controls tension
% % beta: Controls rigidity
% % gamma: Step size
% % kappa: Controls enegry term
% % wl, we, wt: Weights for line, edge and terminal enegy components
% % iterations: No. of iteration for which sn
% xs = col; ys = row;
% P = [xs(:), ys(:)];
% P = MakeContourClockwise2D(P);
% xs = P(:,1)'; ys = P(:,1)';
% alpha = 0.1; beta = 1; gamma = 0.1; kappa = 0.1; wl = 0.1; we = 0.1; wt = 0.1;
% iterations = 50;
% [smth] = interate(im, xs, ys, alpha, beta, gamma, kappa, wl, we, wt, iterations);

%% % implementation 3
% pnts = fliplr(segPoints);
pnts = segPoints;
% pnts = MakeContourClockwise2D(pnts);
alpha = zeros(size(pnts, 1), 1); beta = zeros(size(pnts, 1), 1) + 0.001; 
max_delta_x = 1; resol_x = 1; max_delta_y = 1; resol_y = 1; 
% feat_img is a 2D-Array of the feature responses in the image.  For example it can contain the magnitude of the image gradients
% let's try with just image
feat_img = im;

% now let's try gradient magnitude term
Wline = -1; Wedge = 1;  Wterm = 1; Sigma = 0.01;
Eextern = ExternalForceImage2D(im,Wline, Wedge, Wterm,Sigma);

imshow(im);
hold on
plot(pnts(:,1), pnts(:,2), 'b.');
% recreate the whole segmented area
J = DrawSegmentedArea2D(P,size(im));
LineDisplay.displayMasks(im, J);

[snake_pnts,e] = snake(pnts, alpha, beta, max_delta_x, resol_x, max_delta_y, resol_y, Eextern);

hold on
plot(snake_pnts(:,1), snake_pnts(:,2), 'g.');

