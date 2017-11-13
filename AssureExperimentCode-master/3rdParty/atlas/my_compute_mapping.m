function [mappedX, mapping] = my_compute_mapping(X, type, no_dims)

% This file was originally part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% It has been adapted for our own purposes.
%
% (C) Laurens van der Maaten, Delft University of Technology
% Adaptation: Arturo Mendoza Quispe, Caroline Petitjean, LITIS, Universite
% de Rouen, France - 2015
%
% For the method description, please see (and cite if needed):
% Shape prior based image segmentation using manifold learning
% A. Mendoza Quispe and C. Petitjean, IEEE IPTA, 2015
%
knn = round(size(X,1)/2);
t = 1;
eig = 'Matlab';
sigma = mean(pdist(X));
switch lower(type)
    case 'pca'
        [mappedX, mapping] = compute_mapping(X, type, no_dims);
        mapping.X = X;
    case 'mds'
        [mappedX, mapping] = compute_mapping(X, type, no_dims);
        mapping.X = X;
    case 'isomap'
        [mappedX, mapping] = compute_mapping(X, type, no_dims, knn);
        mapping.X = X;
    case 'lle'
        [mappedX, mapping] = compute_mapping(X, type, no_dims, knn, eig);
        mapping.X = X;
    case {'laplacian', 'laplacianeigenmaps'}
        eig = 'JDQR';
        [mappedX, mapping] = compute_mapping(X, type, no_dims, knn, sigma, eig);
        mapping.X = X;
    case {'dm', 'diffusionmaps'}
        [mappedX, mapping] = my_diffusion_maps(X, no_dims, t, sigma);
        mapping.X = X;
    case 'hessianlle'
        eig = 'JDQR';
        [mappedX, mapping] = compute_mapping(X, type, no_dims, knn, eig);
        mapping.X = X;
    case 'ltsa'
        eig = 'JDQR';
        [mappedX, mapping] = compute_mapping(X, type, no_dims, knn, eig);
        mapping.X = X;
    otherwise; error('Unknown method');
end