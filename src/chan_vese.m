function [seg, phi, hist] = chan_vese(I, phi0, opts)
%CHAN_VESE  Segmentation par level set "sans contours" (Chan & Vese, 2001).
%   [seg, phi, hist] = CHAN_VESE(I, phi0, opts)
%     I    : image en niveaux de gris (double, idealement dans [0,1])
%     phi0 : level set initiale (positive a l'interieur), ou [] pour un
%            cercle centre ; voir init_phi
%     opts : .mu       poids de la longueur du contour       [1]
%            .dt       pas de temps                           [1]
%            .eps      largeur du Heaviside regularise        [1]
%            .nIter    nombre max d'iterations                [800]
%            .reinit   reinitialisation tous les k pas        [10]
%            .tol      arret si moins de cette fraction de pixels change
%                      entre deux reinitialisations (3 fois de suite) [1e-4]
%   seg  : masque logique (phi > 0)
%   phi  : level set finale
%   hist : .phi (cellule de level sets), .energy, .c1, .c2
%
%   Energie :  mu*Longueur(C) + sum_dedans (I - c1)^2 + sum_dehors (I - c2)^2
%   Evolution : dphi/dt = delta_eps(phi) [ mu*div(grad phi/|grad phi|) + F ]
%   avec le terme d'attache aux donnees  -(I - c1)^2 + (I - c2)^2
%                                      = (c1 - c2)(2I - c1 - c2).
%   c1, c2 : moyennes de l'image a l'interieur et a l'exterieur.
%
%   Normalisation : on divise ce terme par (c1 - c2)^2, ce qui donne
%   F = (2I - c1 - c2)/(c1 - c2), egal a +1 sur un pixel typique de
%   l'interieur et -1 a l'exterieur, quel que soit le contraste de l'image ;
%   puis on l'ecrete a [-1, 1] pour qu'un pixel bruite isole ne puisse pas
%   l'emporter sur le terme de longueur. Sans cela, avec une image dans
%   [0,1], le terme de donnees vaut ~0.01 et le contour avance de quelques
%   centiemes de pixel par iteration.
%   Le contour est le niveau zero de phi : il peut se scinder ou fusionner
%   sans traitement particulier (changement de topologie).

if nargin < 2 || isempty(phi0), phi0 = init_phi(size(I), 'cercle'); end
if nargin < 3, opts = struct(); end
mu = getopt(opts, 'mu', 1);
dt = getopt(opts, 'dt', 0.1);
ep = getopt(opts, 'eps', 1);
nIter = getopt(opts, 'nIter', 800);
reinitEvery = getopt(opts, 'reinit', 10);
tol = getopt(opts, 'tol', 1e-4);

phi = phi0;
hist.phi = {phi};
hist.energy = [];
hist.c1 = []; hist.c2 = [];
prevSeg = phi > 0;
stable = 0;
ep_reg=1e-6;
for k = 1:nIter
    Hs = 0.5 * (1 + (2/pi) * atan(phi / ep));
    delta = (ep / pi) ./ (ep^2 + phi.^2);
    c1 = sum(I(:) .* Hs(:)) / max(sum(Hs(:)), eps);
    c2 = sum(I(:) .* (1 - Hs(:))) / max(sum(1 - Hs(:)), eps);
    [gx, gy] = grad_fwd(phi);
    norm= sqrt(gx.^2 + gy.^2 + ep_reg^2);
    nx = gx ./ norm;
    ny = gy ./ norm;
    kappa = div_bwd(nx, ny);
    F = (2*I - c1 - c2) / (c1 - c2 + 1e-12);
    F = min(max(F, -1), 1);
    phi = phi + dt * delta .* (mu * kappa + F);

    % energie normalisee par (c1-c2)^2, comme la force (approximation discrete)

    [gx1, gy1] = gradient(Hs);
    len = sum(sqrt(gx1(:).^2 + gy1(:).^2));
    hist.energy(end+1) = mu * len + (sum((I(:) - c1).^2 .* Hs(:)) + sum((I(:) - c2).^2 .* (1 - Hs(:)))) / (c1 - c2 + 1e-12)^2;
    hist.c1(end+1) = c1; 
    hist.c2(end+1) = c2;
    if mod(k, reinitEvery) == 0
        phi = reinit_sdf(phi, 5);
        hist.phi{end+1} = phi;
        seg = phi > 0;
        change = mean(seg(:) ~= prevSeg(:));
        prevSeg = seg;
        if change < tol
            stable = stable + 1;
            if stable >= 3
                break; 
            end
        else
            stable = 0;
        end
    end
end
seg = phi > 0;
hist.nIter = k;
end

function v = getopt(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
