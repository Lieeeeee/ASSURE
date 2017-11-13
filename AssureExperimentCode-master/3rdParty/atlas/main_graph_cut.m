% (C) Arturo Mendoza Quispe, Caroline Petitjean, LITIS, Universite
% de Rouen, France - 2015
%
% For the method description, please see (and cite if needed):
% Shape prior based image segmentation using manifold learning
% A. Mendoza Quispe and C. Petitjean, IEEE IPTA, 2015
%

function main_graph_cut
global N_LABELS NORM_SIZE;
global BACK_DIM OBJE_DIM OBJE_VAL;
global INI_LABEL;
global FIG_MAIN FIG_TEMP;
N_LABELS = 2;
BACK_DIM = 1;
OBJE_DIM = 2;
OBJE_VAL = 1;
NORM_SIZE = 50;
FIG_MAIN = 1;
FIG_TEMP = 2;

% Parameters
close all;
N_DIMS = 3;
MAX_ITER = 150;
weights = [7/8, 1/8]; % data cost, shape prior cost
max_weights = [5/8, 3/8];
prior_update = 1.05;

% Load data
%%%%%%%%%%%%%%%%%%
% please replace by your own data
%%%%%%%%%%%%%%%%%%%
%distance maps of the learning shape dataset that will be used to construct
%the shape manifold
load(fullfile('Dataset1','data'),'dDistances');
X = dDistances{1}{1}(:,:,1); 

% image to be segmented 
load(fullfile('Dataset1','patient01.mat'),'image');
im = image(:,:,randperm(3,1)); 

% %stack of images of the patient to be segmented 
% required for our initialization - not useful for the general case 
% (please see the initialization function)
load(fullfile('Dataset1','patient01_full.mat'),'image');
ims = image(:,:,80:99); 
clear dDistances image;
%%%%%%%%%%%%%%%%%%

% Initialization of Manifold Learning Procedure
[Y, map] = my_compute_mapping(X, 'DM', N_DIMS);
points = ones(MAX_ITER+1, N_DIMS);
points(1,:) = mean(Y, 1);

% Initialization and placement of shape-prior
norm_shape = reshape(my_pre_image(mean(Y), map), NORM_SIZE, NORM_SIZE);
%The initialization is a in-house process, that is linked to our medical
%data. It might not be useful in your configuration.
[shape_prior, labels, transform] = initialization(im, ims, norm_shape);
INI_LABEL = labels;

% Initialization of Graph Cuts terms
dC = get_data_cost(im, labels);
dC2 = update_data_cost(dC);
sC = get_smoothness_cost();
[vC, hC] = get_smoothness_cost(im);
pC = get_shape_prior_cost(shape_prior, transform);

rep = 0;
weights_ini = weights;
for iter = 1:MAX_ITER
    graphs(iter, weights, im, Y, points, labels, shape_prior, ...
        dC2, pC, norm_shape, transform);
    
    % Segmentation with Shape Prior
    pC = get_shape_prior_cost(shape_prior, transform);
    dC2 = update_data_cost(weights(1)*dC + weights(2)*pC);
    % [vC, hC] = get_smoothness_cost(im, labels);
    gch = GraphCut('open', dC2, sC, vC, hC);
    [gch, labels_new] = GraphCut('swap', gch); % 1 iter
    GraphCut('close', gch);
    
    % Verify conditions
    labels_new = double(labels_new);
    if ~any(labels_new(:)); error('empty'); end
    if any(labels_new(:) - labels(:)); labels = labels_new;
    else rep = rep + 1; end
    if (rep > 10); break; end
    
    % Manifold Transversal
    [shape_prior, transform, points(iter+1,:), norm_shape] = ...
        get_shape_prior(map, labels);
    
    % Update weights
    weights = update_weights(weights, prior_update, max_weights);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = update_weights(w, p, maxw)
w(2) = w(2) * p;
if w(2) > maxw(2); w(2) = maxw(2); end
w(1) = 1 - w(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graphs(iter, weights, im, Y, points, labels, ...
    shape_prior, dC, pC, norm_shape, transformation)
global FIG_MAIN BACK_DIM OBJE_DIM;
figure(FIG_MAIN);
set(gcf,'units','normalized','outerposition', [0 0 1 1])

% Manifold
subplot(2,4,1);
switch size(Y, 2)
    case 2
        plot(Y(:, 1), Y(:, 2), 'b.--', ...
            points(1:iter, 1), points(1:iter, 2), 'rd-', ...
            points(iter, 1), points(iter, 2), 'rs');
    case 3     
        plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'b.--', ...
            points(1:iter, 1), points(1:iter, 2),  points(1:iter, 3), 'rd-', ...
            points(iter, 1), points(iter, 2),  points(iter, 3), 'rs');
    otherwise;
end
grid on; box on; axis equal;
title(sprintf(['Iteration %i\n Current ponderation: E = %.2fE_d + ' ...
    '%.2fE_p + E_s'], iter, weights(1), weights(2)));

% Current Labelling
subplot(2,4,2);
imagesc(im); axis equal; axis off; title('Graph Cuts output');
hold on; h = imagesc(ones(size(im, 1), size(im, 2), 3)); hold off;
set(h, 'AlphaData', 0.5*labels);
subplot(2,8,9);
imagesc(norm_shape); axis equal; axis off;
title('Projected Shape');
subplot(2,8,10);
imagesc(shape_prior); axis equal; axis off;
title('Shape Prior');

% Data Cost
subplot(2,4,3);
imagesc(dC(:, :, OBJE_DIM)); axis equal; axis off;
title('Data Cost Object');
subplot(2,4,7);
imagesc(dC(:, :, BACK_DIM)); axis equal; axis off;
title('Data Cost Background');

% Current Segmentation
subplot(2,4,6);
imagesc(im); axis equal; axis off; title('Segmentation');
hold on; h = imagesc(ones(size(im, 1), size(im, 2), 3)); hold off;
contour = affine_transform(norm_shape <= 0, transformation);
contour = ~bwareaopen(~contour, sum(contour(:))-1);
contour = edge(contour, 'Canny');
set(h, 'AlphaData', contour);
% imagesc(sqrt(vC.^2+hC.^2)); axis equal; axis off; title('smoothness cost');

% Prior Cost
subplot(2,4,4);
imagesc(pC(:, :, OBJE_DIM)); axis equal; axis off;
title('Prior Cost Object');
subplot(2,4,8);
imagesc(pC(:, :, BACK_DIM )); axis equal; axis off;
title('Prior Cost Background');

if (iter < 2); pause; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shape, labels, transformation] = initialization(im, ims, shape)
%The initialization is a in-house process, that is linked to our medical
%cardiac data. It might not be useful in your configuration.
%The aim is just to position the shape prior onto the image.
%You can replace this function with random or manual positioning,
%or any function that fits your data.

global FIG_MAIN NORM_SIZE
done = false;
first = true;
while ~done
    figure(FIG_MAIN);
    set(gcf,'units','normalized','outerposition', [0.05 0.05 0.9 0.9])
    imagesc(im); colormap gray; axis equal; axis off;
    if first
        s = regionprops(find_roi(ims), 'BoundingBox', 'Orientation', 'MinorAxisLength');
        angle = 60 - s.Orientation;
        displacement = s.BoundingBox(1:2)+s.BoundingBox(3:4)*2.1/5;
        d = s.MinorAxisLength*3/5;
        first = false;
    else
        fprintf('Please select the two septum points\n');
        [x, y] = ginput(2);
        angle = atan2(diff(y), diff(x));
        displacement = [mean(x), mean(y)];
        d = sqrt(diff(x)^2 + diff(y)^2)/2;
    end
    T = [1 0 0; 0 1 0; displacement(1) displacement(2) 1];
    R = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
    [m, n] = size(im);
    [xx, yy] = meshgrid(1:n, 1:m);
    origXY = [xx(:), yy(:), ones(numel(xx), 1)];
    [uu, vv] = meshgrid(linspace(-d, d, NORM_SIZE));
    normUV = [uu(:), vv(:), ones(numel(uu), 1)];
    transformation.origXY = origXY;
    transformation.normUV = normUV;
    transformation.matrix = R*T;
    transformation.size = size(im);
    labels = affine_transform(shape <= 0, transformation);
    red = cat(3, ones(size(im)), zeros(size(im, 1), size(im, 2), 2));
    hold on; h = imagesc(red); hold off;
    set(h, 'AlphaData', 0.4*~labels);
    str = input('Is it correct? [Y/N] ', 's');
    figure(FIG_MAIN);
    if strcmp(str, 'y');
        done = true;
    end
end
close(FIG_MAIN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = get_smoothness_cost(varargin)
% sC = get_smoothness_cost()
% [vC, hC] = get_smoothness_cost(im)
% [vC, hC] = get_smoothness_cost(im, labels)
global N_LABELS
switch nargout
    case 1
        sC = 50 * (ones(N_LABELS) - eye(N_LABELS));
        varargout{1} = sC;
    case 2
        im = varargin{1};
        [m, n] = size(im);
        % Smoothness Cost (spatially variant)
        vC = zeros(m, n);
        hC = zeros(m, n);
        varB = 40;
        vf = conv2(fspecial('gauss', [varB, varB], sqrt(varB)), ...
            fspecial('sobel'), 'valid');
        for b=1:size(im,3)
            vC = max(vC, abs(imfilter(im(:,:,b), vf, 'symmetric')));
            hC = max(hC, abs(imfilter(im(:,:,b), vf', 'symmetric')));
        end
        vC = exp(-vC);
        hC = exp(-hC);
        if length(varargin) > 1
            labels = varargin{2};
            vL = double(abs(diff(labels, 1, 1)));
            hL = double(abs(diff(labels, 1, 2)));
            vC = vC.*[vL; zeros(1, n)];
            hC = hC.*[hL, zeros(m, 1)];
        end
        [varargout{1}, varargout{2}] = deal(vC, hC);
    otherwise; error('error');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dC = get_data_cost(im, labels)
global N_LABELS
global BACK_DIM OBJE_DIM OBJE_VAL;
[m, n] = size(im);
dC = zeros(m, n, N_LABELS, 'single');
roi = labels == OBJE_VAL;
se = strel('disk', round(sqrt(sum(roi(:)))*0.25)); % as an square
object = get_seeds(imerode(roi, se));
background = get_seeds(~imdilate(roi, se));
mu(OBJE_DIM) = mean(im(object));
mu(BACK_DIM) = mean(im(background));
sigma(OBJE_DIM) = std(im(object));
sigma(BACK_DIM) = std(im(background));
gaussian = @(x, mu, sigma) exp(-(x-mu).^2/(2*sigma^2));
dC(:, :, OBJE_DIM) = gaussian(im, mu(BACK_DIM), sigma(BACK_DIM));
dC(:, :, OBJE_DIM) = dC(:, :, OBJE_DIM) .* ~object;
dC(:, :, BACK_DIM) = gaussian(im, mu(OBJE_DIM), sigma(OBJE_DIM));
dC(:, :, BACK_DIM) = dC(:, :, BACK_DIM) .* ~background;
dC = update_data_cost(dC);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dC = update_data_cost(dC)
global BACK_DIM OBJE_DIM;
dC(1, :, BACK_DIM) = 0;
dC(end, :, BACK_DIM) = 0;
dC(:, 1, BACK_DIM) = 0;
dC(:, end, BACK_DIM) = 0;
dC(1, :, OBJE_DIM) = 1;
dC(end, :, OBJE_DIM) = 1;
dC(:, 1, OBJE_DIM) = 1;
dC(:, end, OBJE_DIM) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = get_seeds(x)
[m, n] = size(x);
d = 0.05; % density of seeds
s = x(:);
t = length(s(:));
ix = randperm(t, round(t * (1 - d)));
s(ix) = 0;
s = reshape(s, m, n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pC = get_shape_prior_cost(shape, transform)
global N_LABELS;
global BACK_DIM OBJE_DIM;
gamma(OBJE_DIM) = 1.2;
gamma(BACK_DIM) = 25;
maxVal = 1;
pC = zeros(transform.size(1), transform.size(2), N_LABELS);
normalize = @(x) x/max(x(~isinf(x(:)) & ~isnan(x(:))));
T = affine_transform(shape <= 0, transform);
Td = bwdist(T) - bwdist(~T);
M = ones(size(Td));
M2 = ones(size(Td));
M(Td > 0) = exp(-Td(Td > 0)/gamma(OBJE_DIM));
M2(Td > 0) = exp(-Td(Td > 0)/gamma(BACK_DIM));
pC(:, :, BACK_DIM) = normalize(-log(1 - M2)); % Background
pC(:, :, OBJE_DIM) = normalize(-log(M)); % Object
pC(isinf(pC)) = maxVal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shape, transform, point, im] = get_shape_prior(map, labels)
global NORM_SIZE OBJE_VAL
[im, transform] = affine_transform(labels == OBJE_VAL);
point = my_out_of_sample(im(:)', map, 'nystrom');
shape = my_pre_image(point, map);
shape = reshape(shape, NORM_SIZE, NORM_SIZE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = biggest_region(object)
im = bwlabel(object);
h = hist(im(:), max(im(:)) + 1);
[~, ix] = max(h(2:end)); % omit zero
im = im == ix;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function roi = find_roi(ims)
dvals = diff(ims(:,:,2:end),1,3);
roi = log(mean(dvals.^2,3));
roi = roi>mean(roi(:))+std(roi(:));
roi = biggest_region(roi);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = affine_transform(im, transformation)
% [norm, transformation] = normalize(im) ~~ normalization
% t = affine_transform(im, transform) ~~ inverse transform
global NORM_SIZE INI_LABEL
switch nargin
    case 1 % direct transform
        if (nargout ~= 2); error('error'); end
        roi = INI_LABEL;
        se = strel('disk', round(sqrt(sum(roi(:)))*0.25)); % as an square
        roi = imdilate(roi, se);
        im = biggest_region(im & roi);
        [m, n] = size(im);
        stats = regionprops(im, 'Centroid', 'Orientation', 'PixelIdxList');
        [xx, yy] = meshgrid(1:n, 1:m);
        origXY = [xx(:), yy(:), ones(numel(xx), 1)];
        displacement = stats.Centroid;
        angle = deg2rad(stats.Orientation);
        d = max(pdist2(origXY(stats.PixelIdxList, 1:2), displacement));
        d = d * 1.15;
        [uu, vv] = meshgrid(linspace(-d, d, NORM_SIZE));
        normUV = [uu(:) vv(:) ones(numel(uu), 1)];
        T = [1 0 0; 0 1 0; displacement(1) displacement(2) 1];
        R = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];
        normXY = normUV * R * T;
        z = interp2(xx, yy, double(im), ...
            reshape(normXY(:, 1), NORM_SIZE, NORM_SIZE), ...
            reshape(normXY(:, 2), NORM_SIZE, NORM_SIZE), 'nearest');
        z(isnan(z)) = 0;
        norm = bwdist(z) - bwdist(~z);
        transformation.origXY = origXY;
        transformation.normUV = normUV;
        transformation.matrix = R*T;
        transformation.size = size(im);
        varargout{1} = norm;
        varargout{2} = transformation;
    case 2 % inverse transform
        if (nargout ~= 1); error('error'); end
        sz = max(transformation.origXY(:, [2 1]));
        origUV = transformation.origXY / transformation.matrix;
        t = interp2(...
            reshape(transformation.normUV(:, 1), NORM_SIZE, NORM_SIZE), ...
            reshape(transformation.normUV(:, 2), NORM_SIZE, NORM_SIZE), ...
            double(im), ...
            reshape(origUV(:, 1), sz(1), sz(2)), ...
            reshape(origUV(:, 2), sz(1), sz(2)), ...
            'nearest');
        t(isnan(t)) = 0;
        varargout{1} = t;
    otherwise; error('error');
end
end