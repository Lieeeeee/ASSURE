function [ outSeg ] = stapleWrapper( segs )
%STAPLEWRAPPER Summary of this function goes here
%   Detailed explanation goes here


D = zeros(length(segs{1}(:)),length(segs));
for t=1:size(D,2)
   D(:,t) = segs{t}(:); 
end
outSeg = STAPLE(D);
outSeg = reshape(outSeg,size(segs{1},1),size(segs{1},2),size(segs{1},3));
end

