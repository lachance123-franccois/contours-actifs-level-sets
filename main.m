%% Segmentation par contours actifs et level sets
% 1) Snake parametrique de Kass (+ force ballon), schema semi-implicite
% 2) Level set de Chan-Vese (sans contours), avec reinitialisation
% 3) Comparaison avec un seuillage d'Otsu, robustesse au bruit (Dice),
%    changement de topologie, biais d'eclairage, image medicale.
% 4) Etude parametrique : influence de l'initialisation de la level set et
%    de chaque parametre de Chan-Vese et du snake sur le resultat (Dice) et
%    sur le temps de calcul.
% Utilise imgaussfilt (Image Processing Toolbox) pour le lissage de l'image ;
% activecontour sert de reference si elle est presente.
%
% Auteur : Lachance - github.com/lachance123-franccois

clear; close all; clc;
root = fileparts(mfilename('fullpath'));
if ~exist(fullfile(root, 'src'), 'dir'), root = pwd; end   % lance depuis l'editeur
addpath(fullfile(root, 'src'));
resDir = fullfile(root, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
parDir = fullfile(resDir, 'parametres');
if ~exist(parDir, 'dir'), mkdir(parDir); end
rng(3);
sz = 200;

%% 0. Configurations de reference
% Chaque etude parametrique part de ces valeurs et n'en change qu'une a la fois.
% Chan-Vese
refCV.init  = 'sinus';      % type d'initialisation de phi (init_phi)
refCV.sigma =1.5;          % lissage gaussien de l'image avant Chan-Vese (0 = aucun)
refCV.opts  = struct('mu', 0.001, 'dt', 1, 'eps', 1, 'nIter', 4000, 'reinit', 10, 'tol', 1e-4);
% Snake (valeurs par defaut de snake_kass)
refSN = struct('alpha', 0.05, 'beta', 0.1, 'gamma', 1, 'kappa', 2, ...
               'balloon', -0.3, 'sigma', 2, 'nIter', 400, 'nPts', 120);

% Appel de chauffe : le premier appel d'une fonction inclut sa compilation
% (JIT) et fausserait le premier chronometrage.
chan_vese(rand(40), [], struct('nIter', 20));
snake_kass(rand(40), 20 + 10*cos((1:20)'/20*2*pi), 20 + 10*sin((1:20)'/20*2*pi), struct('nIter', 20, 'nPts', 20));

%% 1. Snake vs Chan-Vese vs Otsu sur une image bruitee a deux objets
[I, gt] = make_shapes(sz, 0.15);
%SI = imgaussfilt(I, refCV.sigma); % image bruitee (non filtree)
th = linspace(0, 2*pi, 121)'; th(end) = [];
x0 = sz/2 + 0.45*sz*cos(th); y0 = sz/2 + 0.45*sz*sin(th);
tic; [xs, ys, histS] = snake_kass(I, x0, y0, refSN); tS = toc;
dS = seg_tools('dice', seg_tools('poly2mask', xs, ys, size(I)), gt);
I = imgaussfilt(I, refCV.sigma);                % image filtree, utilisee ensuite
tic; [segCV, phi, histCV] = chan_vese(I, init_phi(size(I), refCV.init), refCV.opts); tC = toc;
dC = seg_tools('dice', seg_tools('bright', segCV, I), gt);
segO = I > seg_tools('otsu', I);
dO = seg_tools('dice', segO, gt);
fprintf('Dice (bruit 0.15) : snake %.3f (%.2f s) | Chan-Vese %.3f (%d it., %.2f s) | Otsu %.3f\n', ...
    dS, tS, dC, histCV.nIter, tC, dO);

fig = figure('Position', [50 50 1300 380]);
subplot(1, 4, 1); imagesc(I); axis image off; colormap gray; hold on;
plot([x0; x0(1)], [y0; y0(1)], 'y--');
title('Image + contour initial');
subplot(1, 4, 2); imagesc(I); axis image off; hold on;
for k = 1:4:numel(histS), c = histS{k}; plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], 'c-'); end
plot([xs; xs(1)], [ys; ys(1)], 'r-', 'LineWidth', 2);
title(sprintf('Snake : Dice %.3f, %.2f s', dS, tS));
subplot(1, 4, 3); imagesc(I); axis image off; hold on;
contour(phi, [0 0], 'r', 'LineWidth', 2);
title(sprintf('Chan-Vese : Dice %.3f, %.2f s', dC, tC));
subplot(1, 4, 4); imagesc(segO); axis image off;
title(sprintf('Otsu : Dice %.3f', dO));
save_fig(fig, fullfile(resDir, 'fig1_comparaison.png'));

%% 2. Changement de topologie : un seul contour initial, deux objets
fig = figure('Position', [50 50 1300 330]);
steps = round(linspace(1, numel(histCV.phi), 5));
for k = 1:5
    subplot(1, 5, k); imagesc(I); axis image off; colormap gray; hold on;
    contour(histCV.phi{steps(k)}, [0 0], 'r', 'LineWidth', 1.5);
    title(sprintf('Iteration %d', (steps(k)-1) * refCV.opts.reinit));
end
sgtitle('Chan-Vese : evolution du contour (niveau zero de \phi)');
save_fig(fig, fullfile(resDir, 'fig2_topologie.png'));

% Animation GIF de l'evolution
gifFile = fullfile(resDir, 'chan_vese_evolution.gif');
fig = figure('Position', [100 100 400 400], 'Color', 'w');
for k = 1:numel(histCV.phi)
    clf; imagesc(I); axis image off; colormap gray; hold on;
    contour(histCV.phi{k}, [0 0], 'r', 'LineWidth', 2);
    title(sprintf('Chan-Vese, iteration %d', (k-1) * refCV.opts.reinit));
    drawnow;
    fr = getframe(fig);
    [A, map] = rgb2ind(fr.cdata, 64);
    if k == 1
        imwrite(A, map, gifFile, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end
close(fig);

%% 3. Energie de Chan-Vese et moyennes c1, c2
fig = figure('Position', [50 50 900 350]);
subplot(1, 2, 1); semilogy(histCV.energy, 'LineWidth', 1.5); grid on;
xlabel('Iteration'); ylabel('Energie / (c_1 - c_2)^2'); title('Energie normalisee');
subplot(1, 2, 2); semilogy([histCV.c1; histCV.c2]', 'LineWidth', 1.5); grid on;
xlabel('Iteration'); legend('c_1 (interieur)', 'c_2 (exterieur)'); title('Moyennes des deux regions');
save_fig(fig, fullfile(resDir, 'fig3_energie.png'));

%% 4. Robustesse au bruit : Dice moyen sur 10 tirages
noise = [0.05 0.1 0.2 0.3 0.4 0.5];
nRep = 10;
D = zeros(numel(noise), 3);
for i = 1:numel(noise)
    for r = 1:nRep
        [In, g] = make_shapes(sz, noise(i));
        sC = seg_tools('bright', chan_vese(In, init_phi(size(In), 'bulles', 20)), In);
        [xx, yy] = snake_kass(In, x0, y0);
        sS = seg_tools('poly2mask', xx, yy, size(In));
        sO = In > seg_tools('otsu', In);
        D(i, :) = D(i, :) + [seg_tools('dice', sS, g), seg_tools('dice', sC, g), seg_tools('dice', sO, g)] / nRep;
    end
    fprintf('bruit %.2f : Dice snake %.3f, Chan-Vese %.3f, Otsu %.3f\n', noise(i), D(i, 1), D(i, 2), D(i, 3));
end
fig = figure('Position', [50 50 700 420]);
plot(noise, D, 'o-', 'LineWidth', 1.5); grid on; ylim([0 1]);
xlabel('Ecart-type du bruit (contraste objet/fond = 0.5)'); ylabel('Dice');
legend('Snake (1 contour pour 2 objets)', 'Chan-Vese', 'Otsu', 'Location', 'southwest');
title(sprintf('Robustesse au bruit (%d tirages par point)', nRep));
save_fig(fig, fullfile(resDir, 'fig4_dice_vs_bruit.png'));

%% 5. Cas d'echec : biais d'eclairage et sensibilite a l'initialisation
[Ib, gb] = make_shapes(sz, 0.05, 0.6);
sOb = Ib > seg_tools('otsu', Ib);
[sCb, phiB] = chan_vese(Ib, init_phi(size(Ib), 'bulles', 20));
sCb = seg_tools('bright', sCb, Ib);
[It, gtt] = make_shapes(sz, 0.15, 0, 'etoile');
[sT1, phiT1] = chan_vese(It, init_phi(size(It), 'cercle', [100 100 90]));
[sT2, phiT2] = chan_vese(It, init_phi(size(It), 'cercle', [30 30 15]), struct('nIter', 150));
fig = figure('Position', [50 50 1300 380]);
subplot(1, 4, 1); imagesc(Ib); axis image off; colormap gray; title('Biais d''eclairage lineaire');
subplot(1, 4, 2); imagesc(sOb); axis image off;
title(sprintf('Otsu : Dice %.2f', seg_tools('dice', sOb, gb)));
subplot(1, 4, 3); imagesc(Ib); axis image off; hold on; contour(phiB, [0 0], 'r', 'LineWidth', 1.5);
title(sprintf('Chan-Vese : Dice %.2f', seg_tools('dice', sCb, gb)));
subplot(1, 4, 4); imagesc(It); axis image off; hold on;
contour(phiT1, [0 0], 'r', 'LineWidth', 1.5); contour(phiT2, [0 0], 'c', 'LineWidth', 1.5);
title(sprintf('Init. large (rouge) %.2f / petite (cyan, 150 it.) %.2f', ...
    seg_tools('dice', seg_tools('bright', sT1, It), gtt), seg_tools('dice', seg_tools('bright', sT2, It), gtt)));
save_fig(fig, fullfile(resDir, 'fig5_limites.png'));

%% 6. Reference : activecontour 
if exist('activecontour', 'file')
    m0 = false(size(I)); m0(20:end-20, 20:end-20) = true;
    sRef = activecontour(I, m0, 500, 'Chan-Vese');
    fprintf('nactivecontour (Chan-Vese, MATLAB) : Dice %.3f\n', seg_tools('dice', sRef, gt));
end

%% 7. Image medicale (optionnel)
% Placez une coupe en niveaux de gris dans data/medical.png et, si vous
% l'avez, la segmentation de reference dans data/medical_gt.png.
mFile = fullfile(root, 'data', 'medical.png');
if exist(mFile, 'file')
    M = double(imread(mFile));
    if size(M, 3) == 3, M = mean(M, 3); end
    M = (M - min(M(:))) / (max(M(:)) - min(M(:)));
    [sM, phiM] = chan_vese(M, init_phi(size(M), 'bulles', 25));
    fig = figure('Position', [50 50 900 420]);
    subplot(1, 2, 1); imagesc(M); axis image off; colormap gray; title('Image medicale');
    subplot(1, 2, 2); imagesc(M); axis image off; hold on; contour(phiM, [0 0], 'r', 'LineWidth', 1.2);
    ttl = 'Chan-Vese';
    gFile = fullfile(root, 'data', 'medical_gt.png');
    if exist(gFile, 'file')
        G = imread(gFile); G = G(:, :, 1) > 127;
        ttl = sprintf('Chan-Vese : Dice %.3f', max(seg_tools('dice', sM, G), seg_tools('dice', ~sM, G)));
    end
    title(ttl);
    save_fig(fig, fullfile(resDir, 'fig6_medical.png'));
else
    fprintf('\n(Ajoutez data/medical.png pour tester sur une image medicale)\n');
end

%% =======================================================================
%  ETUDE PARAMETRIQUE
%  Chan-Vese : image bruitee I0 (bruit 0.15, deux objets).
%  Snake : image "etoile" (un seul objet).
%  Chaque etude part de la configuration de reference (section 0) et ne
%  fait varier qu'un parametre a la fois. Le temps de Chan-Vese inclut le
%  lissage gaussien.
%  =======================================================================
synth = {};                                      % synthese : une ligne par parametre

%% 8. Initialisations de la level set : fonction de niveau et champ
inits = {'cercle',       [];                     % cercle centre (rayon par defaut)
         'cercle',       [sz/2 sz/2 0.45*sz];    % grand cercle englobant
         'cercle',       [0.15*sz 0.15*sz 0.075*sz];   % petit cercle decentre
         'bulles',       20;                     % grille de bulles
         'sinus',        []};                    % damier sinusoidal
noms  = {'cercle', 'grand cercle', 'petit cercle', 'bulles', 'sinus'};
nI = numel(noms);
phi0s = cell(1, nI); phis = cell(1, nI);
dI = zeros(1, nI); tI = zeros(1, nI); kI = zeros(1, nI);
fprintf('\n===== Initialisation de phi (Chan-Vese, autres parametres de reference) =====\n');
for k = 1:nI
    if isempty(inits{k, 2})
        phi0s{k} = init_phi([sz sz], inits{k, 1});
    else
        phi0s{k} = init_phi([sz sz], inits{k, 1}, inits{k, 2});
    end
    tic;
    If = imgaussfilt(I, refCV.sigma);
    [sg, phis{k}, h] = chan_vese(If, phi0s{k}, refCV.opts);
    tI(k) = toc;
    dI(k) = seg_tools('dice', seg_tools('bright', sg, If), gt);
    kI(k) = h.nIter;
    fprintf('%-14s : Dice %.3f | %5d it. | %6.2f s\n', noms{k}, dI(k), kI(k), tI(k));
end
synth(end+1, :) = {'Chan-Vese init', min(dI), max(dI), max(dI) - min(dI), min(tI), max(tI), max(tI) / min(tI)};

[Xg, Yg] = meshgrid(1:sz);
st = 10;                                         % pas du sous-echantillonnage du champ
fig = figure('Position', [20 20 1500 900]);
for k = 1:nI
    % ligne 1 : fonction de niveau initiale phi0 en 3D et plan phi = 0
    subplot(3, nI, k);
    surf(Xg, Yg, phi0s{k}, 'EdgeColor', 'none'); hold on;
    surf(Xg, Yg, zeros(sz), 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    axis tight; view(-35, 35); colormap(gca, parula);
    title(sprintf('\\phi_0 : %s', noms{k}));
    % ligne 2 : lignes de niveau de phi0 et contour initial (niveau zero)
    subplot(3, nI, nI + k);
    imagesc(phi0s{k}); axis image off; colormap(gca, parula); hold on;
    contour(phi0s{k}, 12, 'k');
    contour(phi0s{k}, [0 0], 'r', 'LineWidth', 2);
    title('Lignes de niveau de \phi_0');
    % ligne 3 : phi finale, champ des normales grad(phi)/|grad(phi)|, contour final
    subplot(3, nI, 2*nI + k);
    imagesc(phis{k}); axis image off; colormap(gca, parula); hold on;
    [px, py] = gradient(phis{k});
    nn = hypot(px, py) + eps;
    quiver(Xg(1:st:end, 1:st:end), Yg(1:st:end, 1:st:end), ...
        px(1:st:end, 1:st:end) ./ nn(1:st:end, 1:st:end), ...
        py(1:st:end, 1:st:end) ./ nn(1:st:end, 1:st:end), 0.5, 'w');
    contour(phis{k}, [0 0], 'r', 'LineWidth', 2);
    title(sprintf('\\phi finale : Dice %.3f, %.1f s', dI(k), tI(k)));
end
sgtitle('Influence de l''initialisation : \phi_0 (3D), lignes de niveau, \phi finale et champ \nabla\phi/|\nabla\phi|');
save_fig(fig, fullfile(resDir, 'fig7_initialisations.png'));

%% 9. Effet de la reinitialisation sur la fonction de niveau
% reinit = Inf : jamais de reinitialisation (le test d'arret, fait au moment
% des reinitialisations, est alors desactive). tol = 0 : meme nombre
% d'iterations dans les deux cas.
nR = 300;
oA = refCV.opts; oA.nIter = nR; oA.tol = 0;
oB = oA; oB.reinit = Inf;
If = imgaussfilt(I, refCV.sigma);
p0 = init_phi([sz sz], 'cercle');
tic; [sA, phA] = chan_vese(If, p0, oA); tA = toc;
tic; [sB, phB] = chan_vese(If, p0, oB); tB = toc;
dA = seg_tools('dice', seg_tools('bright', sA, If), gt);
dB = seg_tools('dice', seg_tools('bright', sB, If), gt);
fprintf('\nReinitialisation : avec %.3f (%.2f s) | sans %.3f (%.2f s)\n', dA, tA, dB, tB);
fig = figure('Position', [50 50 1200 650]);
lab = {sprintf('Avec reinitialisation (tous les %d pas)', oA.reinit), 'Sans reinitialisation'};
PH = {phA, phB}; DD = [dA dB]; TT = [tA tB];
for k = 1:2
    [px, py] = grad_fwd(PH{k});
    g = hypot(px, py);
    subplot(2, 3, 3*(k-1) + 1); imagesc(PH{k}); axis image off; colormap(gca, parula); colorbar;
    hold on; contour(PH{k}, [0 0], 'r', 'LineWidth', 1.5);
    title({lab{k}, sprintf('\\phi apres %d it. : Dice %.3f, %.2f s', nR, DD(k), TT(k))});
    subplot(2, 3, 3*(k-1) + 2); imagesc(g, [0 3]); axis image off; colormap(gca, hot); colorbar;
    title('|\nabla\phi| (1 = distance signee)');
    subplot(2, 3, 3*(k-1) + 3); histogram(g(abs(PH{k}) < 10), 50); grid on;
    xlabel('|\nabla\phi| pres du contour (|\phi| < 10)'); ylabel('Nombre de pixels');
end
save_fig(fig, fullfile(resDir, 'fig8_reinitialisation.png'));

%% 10. Chan-Vese : influence de chaque parametre (un seul varie a la fois)
cvSweeps = {
    'sigma',  [0 0.5 1 1.5 2 3],          'Lissage gaussien \sigma (px)';
    'mu',     [0 0.001 0.01 0.1 0.5 1 2], 'Poids de la longueur \mu';
    'dt',     [0.25 0.5 1 2 4],           'Pas de temps dt';
    'eps',    [0.25 0.5 1 2 4],           'Largeur du Heaviside \epsilon';
    'reinit', [2 5 10 20 50],             'Reinitialisation tous les k pas';
    'tol',    [1e-2 1e-3 1e-4 1e-5],      'Tolerance d''arret';
    'nIter',  [50 100 200 500 1000 4000], 'Nombre maximal d''iterations'};
p0 = init_phi([sz sz], refCV.init);
for s = 1:size(cvSweeps, 1)
    name = cvSweeps{s, 1}; vals = cvSweeps{s, 2}; lbl = cvSweeps{s, 3};
    n = numel(vals);
    d = zeros(1, n); t = zeros(1, n); it = zeros(1, n); ph = cell(1, n);
    for k = 1:n
        o = refCV.opts; sig = refCV.sigma;
        if strcmp(name, 'sigma'), sig = vals(k); else, o.(name) = vals(k); end
        tic;
        if sig > 0, If = imgaussfilt(I, sig); else, If = I; end
        [sg, ph{k}, h] = chan_vese(If, p0, o);
        t(k) = toc;
        d(k) = seg_tools('dice', seg_tools('bright', sg, If), gt);
        it(k) = h.nIter;
    end
    % tableau console
    fprintf('\nChan-Vese %s\n   valeur      Dice    iterations   temps (s)\n', name);
    for k = 1:n, fprintf('%9g    %6.3f    %8d    %8.2f\n', vals(k), d(k), it(k), t(k)); end
    synth(end+1, :) = {['Chan-Vese ' name], min(d), max(d), max(d) - min(d), min(t), max(t), max(t) / min(t)}; %#ok<SAGROW>
    % figure : segmentation pour chaque valeur, puis Dice et temps
    fig = figure('Position', [20 60 max(260*n, 900) 640]);
    for k = 1:n
        subplot(2, n, k); imagesc(I); axis image off; colormap(gca, gray); hold on;
        contour(ph{k}, [0 0], 'r', 'LineWidth', 1.5);
        title({sprintf('%s = %g', name, vals(k)), sprintf('Dice %.3f | %d it. | %.2f s', d(k), it(k), t(k))}, 'FontSize', 8);
    end
    subplot(2, 1, 2);
    yyaxis left;  plot(1:n, d, 'o-', 'LineWidth', 1.5); ylabel('Dice'); ylim([0 1]);
    yyaxis right; plot(1:n, t, 's-', 'LineWidth', 1.5); ylabel('Temps (s)');
    set(gca, 'XTick', 1:n, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), vals, 'UniformOutput', false));
    xlim([0.5 n + 0.5]); grid on; xlabel(lbl); title(['Chan-Vese, influence de : ' lbl]);
    save_fig(fig, fullfile(parDir, ['cv_' name '.png'])); close(fig);
end

%% 11. Snake : influence de chaque parametre (image "etoile", un seul objet)
[Ie, ge] = make_shapes(sz, 0.15, 0, 'etoile');
snSweeps = {
    'alpha',   [0.005 0.05 0.2 0.5 1],  'Elasticite \alpha';
    'beta',    [0 0.1 0.5 1 2],         'Rigidite \beta';
    'gamma',   [0.5 1 2 5],             'Pas \gamma (inverse de la viscosite)';
    'kappa',   [0 0.5 2 5 10],          'Poids de la force image \kappa';
    'balloon', [0 -0.1 -0.3 -0.6 -1],   'Force ballon';
    'sigma',   [0.5 1 2 4 6],           'Lissage de la carte de contours \sigma';
    'nIter',   [50 100 200 400 800],    'Nombre d''iterations';
    'nPts',    [30 60 120 240],         'Nombre de points du contour'};
for s = 1:size(snSweeps, 1)
    name = snSweeps{s, 1}; vals = snSweeps{s, 2}; lbl = snSweeps{s, 3};
    n = numel(vals);
    d = zeros(1, n); t = zeros(1, n); xs_ = cell(1, n); ys_ = cell(1, n);
    for k = 1:n
        o = refSN; o.(name) = vals(k);
        tic; [xs_{k}, ys_{k}] = snake_kass(Ie, x0, y0, o); t(k) = toc;
        d(k) = seg_tools('dice', seg_tools('poly2mask', xs_{k}, ys_{k}, size(Ie)), ge);
    end
    it = refSN.nIter * ones(1, n);
    if strcmp(name, 'nIter'), it = vals; end
    fprintf('\nSnake %s\n   valeur      Dice    iterations   temps (s)\n', name);
    for k = 1:n, fprintf('%9g    %6.3f    %8d    %8.2f\n', vals(k), d(k), it(k), t(k)); end
    synth(end+1, :) = {['Snake ' name], min(d), max(d), max(d) - min(d), min(t), max(t), max(t) / min(t)}; %#ok<SAGROW>
    fig = figure('Position', [20 60 max(260*n, 900) 640]);
    for k = 1:n
        subplot(2, n, k); imagesc(Ie); axis image off; colormap(gca, gray); hold on;
        plot([xs_{k}; xs_{k}(1)], [ys_{k}; ys_{k}(1)], 'r-', 'LineWidth', 1.5);
        title({sprintf('%s = %g', name, vals(k)), sprintf('Dice %.3f | %d it. | %.2f s', d(k), it(k), t(k))}, 'FontSize', 8);
    end
    subplot(2, 1, 2);
    yyaxis left;  plot(1:n, d, 'o-', 'LineWidth', 1.5); ylabel('Dice'); ylim([0 1]);
    yyaxis right; plot(1:n, t, 's-', 'LineWidth', 1.5); ylabel('Temps (s)');
    set(gca, 'XTick', 1:n, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), vals, 'UniformOutput', false));
    xlim([0.5 n + 0.5]); grid on; xlabel(lbl); title(['Snake, influence de : ' lbl]);
    save_fig(fig, fullfile(parDir, ['snake_' name '.png'])); close(fig);
end

%% 12. Synthese : quels parametres comptent le plus ?
S = cell2table(synth, 'VariableNames', {'Parametre', 'Dice_min', 'Dice_max', 'Ecart_Dice', ...
    'Temps_min_s', 'Temps_max_s', 'Rapport_temps'});
S = sortrows(S, 'Ecart_Dice', 'descend');
fprintf('\n===== Synthese : sensibilite du Dice et du temps a chaque parametre =====\n');
disp(S);
writetable(S, fullfile(parDir, 'synthese_parametres.csv'));
fprintf('\nFigures et GIF enregistres dans %s\n', resDir);
fprintf('Etude parametrique dans %s\n', parDir);