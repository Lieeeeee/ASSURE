function [ nPoints ] = getNPointsBetweenP1P2(x1,y1,x2,y2)
%GETNPOINTSBETWEENP1P2 Summary of this function goes here
%   Detailed explanation goes here

nPoints = ceil(sqrt((x2 - x1).^2 + (y2 - y1).^2 + 0.0000001));

end

