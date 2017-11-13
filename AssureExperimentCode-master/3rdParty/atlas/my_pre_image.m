function recX = my_pre_image(mappedX, mapping)
%MY_PRE_IMAGE Reconstructs data from low-dimensional data representation
%
%   recX = my_pre_image(mappedX, mapping)
%
% The function reconstructs the data in its original space from its
% low-dimensional representation stored in mappedX, using Delaunay triangulation. The function requires the
% original mapping as input. The reconstructed data is stored in recX. 
%
% For the method description, please see (and cite if needed):
% Shape prior based image segmentation using manifold learning
% A. Mendoza Quispe and C. Petitjean, IEEE IPTA, 2015
%
% (C) Arturo Mendoza Quispe, Caroline Petitjean, LITIS, Universite
% de Rouen, France - 2015
%


if ~strcmp(mapping.name, 'DM')
    recX = reconstruct_data(mappedX, mapping);
    X = mapping.X;
    n = min(X, [], 1);
    m = max(X, [], 1);
    recX(recX > m) = m(recX > m);
    recX(recX < n) = m(recX < n);
else
    alpha = 0.5; % gradient descent
    normal = @(x) x/sum(x);
    d = size(mapping.vec, 2);
    tri = delaunayTriangulation(mapping.vec);
    recX = zeros(size(mappedX, 1), size(mapping.X, 2));
    neighbors = zeros(size(mappedX, 1), d+1);
    X = mapping.X;
    maxIter = 1e2;
    delta = 1e-3;
    for i = 1:size(mappedX, 1)
        t_point = mappedX(i, :);
        [ti, theta] = pointLocation(tri, t_point);
        if isnan(ti)
            t_point = tri.Points(nearestNeighbor(tri, t_point), :);
            [ti, theta] = pointLocation(tri, t_point);
        end
        t_neighbors = tri(ti, :)';
        theta = theta(1:d+1)';
        for iter = 1:maxIter
            theta_old = theta;
            t_pre_image = sum(bsxfun(@times, X(t_neighbors,:), theta_old), 1);
            dist = sum(bsxfun(@minus, X(t_neighbors,:) <= 0, ...
                t_pre_image <= 0).^2, 2);
            % J(i) = sum(theta.*dist);
            theta = normal(alpha * normal(dist) + (1 - alpha) * theta_old);
            if (sum(abs(theta./theta_old - 1)) < delta); break; end;
        end
        recX(i,:) = t_pre_image + correction;
        neighbors(i,:) = t_neighbors(:)';
    end
end