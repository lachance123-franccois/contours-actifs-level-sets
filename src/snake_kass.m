function [x, y, hist] = snake_kass(I, x0, y0, opts)
%SNAKE_KASS  Contour actif parametrique (Kass, Witkin & Terzopoulos, 1988).
%   [x, y, hist] = SNAKE_KASS(I, x0, y0, opts)
%     I      : image (double)
%     x0, y0 : points du contour initial ferme (pixels)
%     opts   : .alpha   elasticite (tension)             [0.05]
%              .beta    rigidite (courbure)               [0.1]
%              .gamma   pas (inverse de la viscosite)     [1]
%              .kappa   poids de la force image           [2]
%              .balloon force normale (<0 : retrecit)     [-0.3]
%              .sigma   lissage avant la carte de contours [2]
%              .nIter   iterations                        [400]
%              .nPts    nombre de points du contour       [120]
%   x, y : contour final ; hist : contours intermediaires (cellules)
%
%   Energie : sum alpha|v'|^2 + beta|v''|^2 - kappa |grad(G_sigma * I)|^2
%   Schema semi-implicite : (A + gamma Id) v_{t+1} = gamma v_t + F_ext(v_t)
%   ou A est la matrice pentadiagonale circulante de l'energie interne.
%   La force "ballon" (Cohen, 1991) pousse le contour le long de sa normale,
%   ce qui permet de partir loin de l'objet.

if nargin < 4, opts = struct(); end
alpha = getopt(opts, 'alpha', 0.05);
beta = getopt(opts, 'beta', 0.1);
gamma = getopt(opts, 'gamma', 1);
kappa = getopt(opts, 'kappa', 2);
balloon = getopt(opts, 'balloon', -0.3);
sigma = getopt(opts, 'sigma', 2);
nIter = getopt(opts, 'nIter', 400);
nPts = getopt(opts, 'nPts', 120);

% Carte de contours et force externe normalisee
h = ceil(3*sigma); t = -h:h;
g = exp(-t.^2 / (2*sigma^2)); g = g / sum(g);
Is = conv2(g, g, padarray_rep(I, h), 'valid');
[gx, gy] = gradient(Is);
Eedge = gx.^2 + gy.^2;
Eedge = Eedge / max(Eedge(:));
[fx, fy] = gradient(Eedge);
m = max(hypot(fx(:), fy(:)));
fx = fx / m; fy = fy / m;

% Matrice interne (circulante) et son inverse
[x, y] = resample_closed(x0(:), y0(:), nPts);
n = numel(x);
a = beta; b = -alpha - 4*beta; c = 2*alpha + 6*beta;
A = c * eye(n) + b * (circshift(eye(n), 1) + circshift(eye(n), -1)) ...
    + a * (circshift(eye(n), 2) + circshift(eye(n), -2));
Minv = inv(A + gamma * eye(n));

hist = {[x y]};
[H, W] = size(I);
for k = 1:nIter
    Fx = interp2(fx, x, y, 'linear', 0);
    Fy = interp2(fy, x, y, 'linear', 0);
    % normale (sens antihoraire -> normale exterieure)
    tx = circshift(x, -1) - circshift(x, 1);
    ty = circshift(y, -1) - circshift(y, 1);
    nt = hypot(tx, ty) + eps;
    nx = ty ./ nt; ny = -tx ./ nt;
    if signed_area(x, y) < 0, nx = -nx; ny = -ny; end
    x = Minv * (gamma * x + kappa * Fx + balloon * nx);
    y = Minv * (gamma * y + kappa * Fy + balloon * ny);
    x = min(max(x, 1), W); y = min(max(y, 1), H);
    if mod(k, 10) == 0
        [x, y] = resample_closed(x, y, nPts);   % repartition uniforme des points
        hist{end+1} = [x y]; 
    end
end
end

function J = padarray_rep(I, h)
[H, W] = size(I);
J = I(min(max((1-h):(H+h), 1), H), min(max((1-h):(W+h), 1), W));
end

function A = signed_area(x, y)
A = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y);
end

function [xr, yr] = resample_closed(x, y, n)
xc = [x; x(1)]; yc = [y; y(1)];
s = [0; cumsum(hypot(diff(xc), diff(yc)))];
[s, iu] = unique(s);
xc = xc(iu); yc = yc(iu);
si = linspace(0, s(end), n + 1)'; si(end) = [];
xr = interp1(s, xc, si, 'linear');
yr = interp1(s, yc, si, 'linear');
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
