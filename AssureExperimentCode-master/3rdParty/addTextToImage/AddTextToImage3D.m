function outImage = AddTextToImage3D(Image,String,Position,Color,Font,FontSize)

%rgbImage
outImage = zeros(size(Image));

%case RGB+depth
if(size(outImage,4))==3
    for z=1:size(Image,3)
        inputIm = squeeze(Image(:,:,z,:));
        
        if ~exist('Color','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position);
        elseif ~exist('Font','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color);
        elseif ~exist('FontSize','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color,Font);
        else
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color,Font,FontSize);
        end
        
        outImage(:,:,z,:) =  currentOut;
    end
    
else
    %case grayscale+depth
    for z=1:size(Image,3)
        inputIm = squeeze(Image(:,:,z));
        if ~exist('Color','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position);
        elseif ~exist('Font','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color);
        elseif ~exist('FontSize','var')
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color,Font);
        else
            currentOut = AddTextToImage(squeeze(inputIm),String,Position,Color,Font,FontSize);
        end
        outImage(:,:,z) =  currentOut;
    end
    
end

end