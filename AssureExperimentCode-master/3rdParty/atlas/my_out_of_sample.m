function t_point = my_out_of_sample(point, mapping, method)
%MY_OUT_OF_SAMPLE Performs out-of-sample extension of new datapoints
%
%   t_points = my_out_of_sample(points, mapping, method)
%
% Two methods are implemented to perform out-of-sample extension of the new datapoints, the Nystrom extension
% and the sparse representation. 
%
% For the method description, please see (and cite if needed):
% Shape prior based image segmentation using manifold learning
% A. Mendoza Quispe and C. Petitjean, IEEE IPTA, 2015
%
%

if ~strcmp(mapping.name, 'DM')
    t_point = out_of_sample(point, mapping);
else
    if (nargin < 3); method = 'nystrom'; end
    t_point = zeros(size(point, 1),numel(mapping.val));
    switch lower(method)
        case 'nystrom'
            for i = 1:size(point, 1)
                pX = double(point(i,:));
                K = exp((-pdist2(mapping.X,pX).^2)/(2*mapping.sigma^2));
                p = sum(K, 1)';
                K = K ./ ((p * p') .^ mapping.t);
                p = sqrt(sum(K, 1))';
                K = K ./ (p * p');
                P = repmat(K, [1 size(mapping.vec, 2)]);
                eigfun = sum(mapping.vec .* P, 1);
                t_point(i,:) = eigfun./mapping.val' + correction;
            end
        case 'sparse'
            [N, D] = size(mapping.X);
            XI = [mapping.X', eye(D)];
            a0 = zeros(N + D, 1);
            a0(randperm(N,4)) = 1/4;
            for i = 1:size(point, 1)
                b = double(point(i,:))';
                a = l1eq_pd(a0, XI, [], b, 1e-3, 5);
                w = a(1:N);
                % e = a(N+1:end);
                t_point(i,:) = sum(bsxfun(@times,mapping.vec,w),1)/sum(w);
            end
        otherwise; error('error');
    end
end