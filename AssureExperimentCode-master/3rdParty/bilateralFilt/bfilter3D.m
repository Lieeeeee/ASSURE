function [ B ] = bfilter3D(A,w,sigma)
%BFILTER3D Summary of this function goes here
%   Detailed explanation goes here
if ~exist('sigma','var')
    sigma = [];
end
if ~exist('w','var')
    w = [];
end
B = zeros(size(A));

for z=1:size(A,3)
   out = bfilter2(A(:,:,z),w,sigma); 
   B(:,:,z) = out;
end

end

