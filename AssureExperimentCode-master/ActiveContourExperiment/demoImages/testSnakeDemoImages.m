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
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg, OptionsIn, OptionsOut);
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
varMask = VariabilityEstimator.evaluate3DVarMaskActiveContours(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);
im0 = LineDisplay.getCroppedFrameFromFigure();
LineDisplay.displayVariabilityFromMask(I, seg, varMask, false);
im1 = LineDisplay.getCroppedFrameFromFigure();
combIm = [im0,im1];
imshow(combIm);

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
