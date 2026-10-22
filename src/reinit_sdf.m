function phi = reinit_sdf(phi, nIter)

if nargin < 2, nIter = 5; end
dt = 0.1;
S = phi ./ sqrt(phi.^2 + 1);
for k = 1:nIter
    a = phi - phi(:, [1 1:end-1]);       % D-x
    b = phi(:, [2:end end]) - phi;       % D+x
    c = phi - phi([1 1:end-1], :);       % D-y
    d = phi([2:end end], :) - phi;       % D+y
    ap = max(a, 0); am = min(a, 0);
    bp = max(b, 0); bm = min(b, 0);
    cp = max(c, 0); cm = min(c, 0);
    dp = max(d, 0); dm = min(d, 0);

    gPos = sqrt(max(ap.^2, bm.^2) + max(cp.^2, dm.^2));
    gNeg = sqrt(max(am.^2, bp.^2) + max(cm.^2, dp.^2));
    G = (S > 0) .* gPos + (S < 0) .* gNeg + (S == 0);
    phi = phi + dt * S .* (1 - G);
end
end
