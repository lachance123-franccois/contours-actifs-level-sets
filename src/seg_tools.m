function varargout = seg_tools(action, varargin)
%SEG_TOOLS  Petits outils de segmentation (sans toolbox).
%   d   = SEG_TOOLS('dice', A, B)        coefficient de Dice 2|A∩B|/(|A|+|B|)
%   t   = SEG_TOOLS('otsu', I)           seuil d'Otsu (variance inter-classes max)
%   M   = SEG_TOOLS('poly2mask', x, y, sz) masque d'un polygone ferme
%   seg = SEG_TOOLS('bright', seg, I)    choisit la phase la plus claire comme objet

switch action
    case 'dice'
        A = varargin{1}; B = varargin{2};
        varargout{1} = 2 * nnz(A & B) / (nnz(A) + nnz(B));
    case 'otsu'
        I = varargin{1};
        edges = linspace(min(I(:)), max(I(:)), 257);
        p = histcounts(I(:), edges);
        p = p / sum(p);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        w0 = cumsum(p);
        mu0 = cumsum(p .* centers);
        muT = mu0(end);
        sb = (muT * w0 - mu0).^2 ./ (w0 .* (1 - w0) + eps);
        [~, k] = max(sb);
        varargout{1} = centers(k);
    case 'poly2mask'
        x = varargin{1}; y = varargin{2}; sz = varargin{3};
        [X, Y] = meshgrid(1:sz(2), 1:sz(1));
        varargout{1} = inpolygon(X, Y, x, y);
    case 'bright'
        seg = varargin{1}; I = varargin{2};
        if mean(I(seg)) < mean(I(~seg)), seg = ~seg; end
        varargout{1} = seg;
    otherwise
        error('seg_tools:action', 'Action inconnue : %s', action);
end
end
