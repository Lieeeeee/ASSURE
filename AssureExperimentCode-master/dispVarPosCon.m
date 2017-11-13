function dispVarPosCon( imgStrct )
%DISPVARPOSCON Summary of this function goes here
%   Detailed explanation goes here

[union,intersection,uncertaintyMask] = Utils.calcUnionIntersection(imgStrct.masks);

LineDisplay.displayMasks(imgStrct.img,{intersection});
im0 = LineDisplay.getCroppedFrameFromFigure();
close all;

LineDisplay.displayMasks(imgStrct.img,{zeros(size(union)),union});
im1 = LineDisplay.getCroppedFrameFromFigure();
close all;

LineDisplay.displayMasks(imgStrct.img,{zeros(size(union)),zeros(size(union)),uncertaintyMask});
im2 = LineDisplay.getCroppedFrameFromFigure();
close all;

LineDisplay.displayMasks(imgStrct.img,{zeros(size(union))});
im3 =  LineDisplay.getCroppedFrameFromFigure();
close all;

LineDisplay.displaySegsOverlay(imgStrct.img,imgStrct.masks,[],false);
im4 =  LineDisplay.getCroppedFrameFromFigure();
close all;

combIm = [im3,im4,im0,im1,im2];
widths = [size(im3,2),size(im4,2),size(im0,2),size(im1,2),size(im2,2)];
imshow(combIm);

[h] = size(im4,1);
for z=1:length(widths)
    text( sum(widths(1:z))-widths(z)/2,h+20,['(' num2str(z) ')'],'FontSize',12,'FontWeight','bold','Color','black');
end
end

