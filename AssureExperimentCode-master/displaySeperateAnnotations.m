function im = displaySeperateAnnotations( imgCell, slice, nAnnotations )
%DISPLAYSEPERATEANNOTATIONS Summary of this function goes here
%   Detailed explanation goes here

for z=1:length(imgCell.masks)
    imgCell.masks{z} =  imgCell.masks{z}(:,:,slice);
end

mat2show = [];
for j=1:nAnnotations
    currMasks =  cell(size(imgCell.masks));
    for z=1:length(currMasks)
       if z~=j
           currMasks{z} = zeros(size(imgCell.masks{1}));
       else
           currMasks{z} = imgCell.masks{j};
       end
    end
    LineDisplay.displaySegsOverlay(imgCell.img(:,:,slice),currMasks);
    legendMat = LineDisplay.getCroppedFrameFromFigure();
    close all;
    mat2show = [mat2show,legendMat];
end
LineDisplay.displaySegsOverlay(imgCell.img(:,:,slice), imgCell.masks(1:nAnnotations));
legendMat = LineDisplay.getCroppedFrameFromFigure();
close all;
mat2show = [mat2show,legendMat];
LineDisplay.displaySegsOverlay(imgCell.img(:,:,slice), {zeros(size(imgCell.masks{1}))});
legendMat = LineDisplay.getCroppedFrameFromFigure();
mat2show = [legendMat,mat2show];
imshow(mat2show);


im = mat2show;
end

