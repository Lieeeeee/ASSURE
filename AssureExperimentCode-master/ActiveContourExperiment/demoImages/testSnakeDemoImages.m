% testing quality of snakes variability estimation using synthetic images
[OptionsIn, OptionsOut] = SnakeOptions.getLungOptions(1);

%% dark circle with aligned segmentation
[I, seg] = DemoImagesCl.getCircleAndSeg();
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im0 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);
im1 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1];
imshow(combIm);

%% dark circle with missaligned segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getDarkCircleImg());
seg = DemoImagesCl.getCircleSeg(15, [50, 55]);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im0 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);
im1 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1];
imshow(combIm);

%% light circle with aligned segmentation
[I, seg] = DemoImagesCl.getCircleAndSeg(0);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);

%% light circle with missaligned segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getLightCircleImg());
seg = DemoImagesCl.getCircleSeg(15, [50, 55]);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);

%% %%%%%%%%%%%%% ACTIVE CONTOURS %%%%%%%%%%%%%
% [OptionsIn, OptionsOut] = ActiveContourOptions.getLungOptions();
outDir = '/cs/casmip/clara.herscu/git/thesis/figs/activeContour_feb18/demo/';

params = [10, 0.01, 0, 10, -0.01, 0]; 
iterIn = params(1); cIn = params(2); sIn = params(3); 
OptionsIn = ActiveContourOptions.getSpecifiedOptions(iterIn, cIn, sIn);

iterOut = params(4); cOut = params(5); sOut = params(6); 
OptionsOut = ActiveContourOptions.getSpecifiedOptions(iterOut, cOut, sOut);

%% dark circle with aligned segmentation
[I, seg] = DemoImagesCl.getCircleAndSeg();
LineDisplay.displayMasks(I, seg);
im0 = LineDisplay.getCroppedFrameFromFigure();
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im1 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);
im2 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1,im2];
close all;
figure, imshow(combIm);

%% dark circle with missaligned segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getDarkCircleImg());
seg_moved = DemoImagesCl.getCircleSeg(15, [50, 55]);
LineDisplay.displayMasks(I, seg_moved);
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg_moved, OptionsIn, OptionsOut);
im0 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im1 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg_moved, varMask, false);
im2 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1,im2];
close all;
figure, imshow(combIm);

%% light circle with aligned segmentation
[I, seg] = DemoImagesCl.getCircleAndSeg(0);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);


%% light circle with missaligned segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getLightCircleImg());
seg = DemoImagesCl.getCircleSeg(15, [50, 55]);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);


%%  dark circle with ellipse segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getDarkCircleImg());
hor_radius = 10; ver_radius = 30; location = [50 50];
seg_ellipse = DemoImagesCl.getEllipseSeg(hor_radius, ver_radius, location);
LineDisplay.displayMasks(I, seg_ellipse);
im0 = LineDisplay.getCroppedFrameFromFigure();
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg_ellipse, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im1 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg_ellipse, varMask, false);
im2 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1,im2];
close all;
figure, imshow(combIm);


%% dark circle with big segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getDarkCircleImg());
seg_big = DemoImagesCl.getCircleSeg(25, [50, 50]);
LineDisplay.displayMasks(I, seg_big);
im0 = LineDisplay.getCroppedFrameFromFigure();
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg_big, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im1 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg_big, varMask, false);
im2 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1,im2];
close all;
figure, imshow(combIm);
