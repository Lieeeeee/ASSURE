function [rgb]=gray2rgb(Image)
%Gives a grayscale image an extra dimension
%in order to use color within it
Image = mat2gray(Image);
[m n]=size(Image);
rgb=zeros(m,n,3);
rgb(:,:,1)=Image;
rgb(:,:,2)=rgb(:,:,1);
rgb(:,:,3)=rgb(:,:,1);
end
