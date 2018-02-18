% testing quality of snakes variability estimation using synthetic images
[OptionsIn, OptionsOut] = SnakeOptions.getLungOptions(1);

%% dark circle with aligned segmentation
[I, seg] = DemoImagesCl.getCircleAndSeg();
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);

%% dark circle with missaligned segmentation
I = DemoImagesCl.addBlurToImg(DemoImagesCl.getDarkCircleImg());
seg = DemoImagesCl.getCircleSeg(15, [50, 55]);
LineDisplay.displayMasks(I, seg);
varMask = VariabilityEstimator.evaluate3DVarMaskSnakes(I, seg, OptionsIn, OptionsOut);
LineDisplay.displayVariabilityWithoutSeg(I, varMask);

