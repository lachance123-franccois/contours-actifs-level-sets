function phi = init_phi(sz, type, param)
%INIT_PHI  Level sets initiales (distance signee, positive a l'interieur).
%   phi = INIT_PHI(sz, 'cercle', [cx cy r])  cercle (coordonnees en pixels)
%                                           defaut : centre, rayon 40 % du cote
%   phi = INIT_PHI(sz, 'bulles', pas)        grille de petites bulles (insensible
%                                           a l'initialisation, defaut pas = 20 px)
H = sz(1); W = sz(2);
[X, Y] = meshgrid(1:W, 1:H);
switch type
    case 'cercle'
        if nargin < 3 || isempty(param), param = [W/2 H/2 0.4*min(H, W)]; end
        phi = param(3) - hypot(X - param(1), Y - param(2));
    case 'bulles'
        if nargin < 3 || isempty(param), param = 20; end
        r = param / 3;
        dx = mod(X, param) - param/2;
        dy = mod(Y, param) - param/2;
        phi = r - hypot(dx, dy);
    case 'sinus'
        if nargin < 3 || isempty(param)
            L = W; 
        else
            L = param(1);
        end
        phi = sin((pi * X) / L) .* sin((pi * Y) / L);

    otherwise
        error('init_phi:type', 'Type inconnu : %s', type);
end
end
