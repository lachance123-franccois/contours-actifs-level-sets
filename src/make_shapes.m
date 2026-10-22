function [I, gt] = make_shapes(sz, noiseStd, inhomog, layout)
%MAKE_SHAPES  Image synthetique d'objets sur fond, avec verite terrain.
%   [I, gt] = MAKE_SHAPES(sz, noiseStd, inhomog, layout)
%     sz       : taille (pixels)                                    [200]
%     noiseStd : ecart-type du bruit gaussien (contraste objet/fond = 0.5) [0.1]
%     inhomog  : amplitude d'un biais d'eclairage lineaire           [0]
%     layout   : 'deux' (disque + carre arrondi) ou 'etoile'        ['deux']
%   I  : image double
%   gt : masque logique des objets

if nargin < 1 || isempty(sz), sz = 200; end
if nargin < 2 || isempty(noiseStd), noiseStd = 0.1; end
if nargin < 3 || isempty(inhomog), inhomog = 0; end
if nargin < 4 || isempty(layout), layout = 'deux'; end

[X, Y] = meshgrid(linspace(0, 1, sz));
switch layout
    case 'deux'
        gt = hypot(X - 0.30, Y - 0.35) < 0.17;
        sq = max(abs(X - 0.68), abs(Y - 0.62));
        gt = gt | (sq < 0.17 & hypot(max(abs(X - 0.68) - 0.10, 0), max(abs(Y - 0.62) - 0.10, 0)) < 0.07);
    case 'etoile'
        th = atan2(Y - 0.5, X - 0.5);
        rr = hypot(X - 0.5, Y - 0.5);
        gt = rr < 0.25 + 0.08 * cos(5 * th);
    otherwise
        error('make_shapes:layout', 'Disposition inconnue : %s', layout);
end
I = 0.25 + 0.5 * gt + inhomog * (X - 0.5) + noiseStd * randn(sz);
end
