function [ K, actual_nPoints ] = SampleContourPoints2D( P,nPoints )
%SAMPLECONTOURPOINTS2D this function samples nPoints on the contour.
%
% K, actual_nPoints =InterpolateContourPoints(P,nPoints)
%
% input,
%  P : Inpute Contour, size N x 2  (with N>=4)
%  nPoints : Number of Contour points as output
% 
% output,
%  K : Uniform sampled Contour points, size actual_nPoints x 2
%  actual_nPoints: number of contour points in the result

stepSize = ceil(size(P, 1) / nPoints);
actual_nPoints = size(P, 1) / stepSize;

K = P(1:stepSize:end, :);

end

