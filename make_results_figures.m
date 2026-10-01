% version v10, 26/09/2026 : F remplace par -log10(p_FDR), critere reel de la
%                          porte 1 avec eta2 ; ordre des panneaux = ordre des portes
% version v9, 20/09/2026 : une barre par panneau, cases NaN transparentes,
%                          precision relative a la place de C
% version v8, 19/09/2026 : matrices sans valeurs, cases hautes, texte a 11 pt
% version v6, 19/09/2026 : matrices plus hautes, panneaux resserres
% version v5, 19/09/2026 : barre de couleur partagee et etiquettes de descripteurs
%                          sur le seul panneau de gauche dans les matrices
% version v4, 19/09/2026 : legende partagee sous les courbes, graduations impaires,
%                          cases proportionnees dans les matrices
% version v3, 31/07/2026 : panel a six descripteurs, D4 correle mais peu fiable
%MAKE_RESULTS_FIGURES  Figures de resultats du volet 1.
%   Produit quatre fichiers :
%     Fig_matrices_discrim.png   -log10(p_FDR), eta2, ICC  (quantites calculees sur X)
%     Fig_matrices_function.png  R2, Spearman, r           (quantites faisant intervenir Y)
%     Fig_curves_discrim.png     memes quantites en courbes
%     Fig_curves_function.png    idem
%
%   POURQUOI p_FDR ET NON F. Avec le branchement, la p-valeur ne vient pas
%   toujours du F standard (Welch ou Kruskal-Wallis selon les hypotheses).
%   p_FDR est la seule quantite commune a toutes les paires, et c'est elle qui
%   entre dans la porte 1 avec eta2. Elle est tracee en -log10, plafonnee a
%   PLAFOND pour que les valeurs numeriquement nulles n'ecrasent pas le reste.
%   La ligne de seuil est a -log10(P.q), soit 1.30 pour q = 0.05.
%
%   Usage :  make_results_figures   (lit ws_run_volet1.mat)
S1 = load('ws_run_volet1.mat','out');
out = S1.out;
P = paper_params();
names = {'D1','D2','D3','D4','D5'};
S = P.S; K = P.K;

% ===================== STYLE =================================================
St.fs_ax    = 11;      % graduations
St.fs_lab   = 12;      % titres d'axes
St.fs_title = 13;      % titres de panneaux
St.fs_matrix= 11;      % texte des matrices, a la taille finale
St.fs_leg   = 10;
St.lw       = 1.8;     % epaisseur des courbes
St.ms       = 6;       % taille des marqueurs
St.dpi      = 400;
St.plafond  = 15;      % plafond d'affichage de -log10(p_FDR)
% =============================================================================

% groupe 1 : quantites des portes 1 et 2, calculees sur X seul, dans l'ordre
% des portes. La 3e colonne 'mlog10' signale la transformation -log10.
G1 = { out.p_fdr, '-log_{10} p_{FDR}', 'mlog10', -log10(P.q) ;
       out.eta2,  '\eta^2',            [0 1],    P.eta2_thr ;
       out.ICC,   'ICC',               [0 1],    P.icc_thr };

% groupe 2 : quantites faisant intervenir la reponse fonctionnelle
G2 = { out.R2m,   'R^2 marginal',         [0 1],   P.r2_thr ;
       out.rho,   'Spearman \rho',         [-1 1],  [] ;
       out.b1_rp, 'Relative precision r',  [0 0.2], [] };
% La coherence de signe C vaut 1 pour toute paire retenue et n'apprend rien sur
% une carte ; la precision relative de la pente la remplace. Elle se lit a
% l'inverse des deux autres : plus basse veut dire mieux.

drawMatrices(G1, 'Fig_matrices_discrim.png',  names, S, K, St);
drawMatrices(G2, 'Fig_matrices_function.png', names, S, K, St);
drawCurves(G1,  'Fig_curves_discrim.png',  names, S, K, St);
drawCurves(G2,  'Fig_curves_function.png', names, S, K, St);
fprintf('Quatre figures sauvees.\n');

% -----------------------------------------------------------------------------
function Y = transformer(D, lim, St)
% -log10 plafonne pour les p-valeurs, identite sinon. Les NaN restent NaN.
if ischar(lim) && strcmp(lim,'mlog10')
    Y = min(-log10(max(D, 10^-St.plafond)), St.plafond);
else
    Y = D;
end
end

% -----------------------------------------------------------------------------
function drawMatrices(G, fname, names, S, K, St)
% Cartes descripteur-echelle, une par quantite, chacune avec sa propre barre de
% couleur. Les paires sans valeur (NaN) sont transparentes sur fond blanc.
n = size(G,1);
fig = figure;
t = tiledlayout(fig, 1, n, 'TileSpacing','compact', 'Padding','compact');
for q = 1:n
    D = G{q,1}; lim = G{q,3};
    ax = nexttile(t);
    Y = transformer(D, lim, St);
    if ischar(lim)
        clims = [0 St.plafond];
    else
        clims = lim;
    end
    hI = imagesc(ax, Y', clims);
    clim(ax, clims); colorbar(ax);
    set(hI, 'AlphaData', ~isnan(Y'));
    set(ax, 'Color', [1 1 1]);
    colormap(ax, parula(256));
    set(ax,'XTick',1:2:S,'XTickLabelRotation',0);
    pbaspect(ax, [S 1.9*K 1]);
    xlabel(ax,'Scale level');
    if q == 1
        set(ax,'YTick',1:K,'YTickLabel',names);
        ylabel(ax,'Descriptor');
    else
        set(ax,'YTick',1:K,'YTickLabel',[]);
    end
    title(ax, G{q,2}, 'FontWeight','bold');
end
paper_figstyle(fig, 15, 8);
bumpFonts(fig, St.fs_matrix);
exportgraphics(fig, fname, 'Resolution', St.dpi); close(fig);
end

% -----------------------------------------------------------------------------
function bumpFonts(fig, fs)
% paper_figstyle harmonise toutes les figures a 10 points. Les matrices portent
% peu de texte et supportent un point de plus, applique ici seulement.
for ax = findall(fig,'Type','axes').'
    set(ax,'FontSize',fs);
    set(get(ax,'XLabel'),'FontSize',fs);
    set(get(ax,'YLabel'),'FontSize',fs);
    set(get(ax,'Title'), 'FontSize',fs,'FontWeight','bold');
end
for cb = findall(fig,'Type','colorbar').'
    set(cb,'FontSize',fs);
end
end

% -----------------------------------------------------------------------------
function drawCurves(G, fname, names, S, K, St)
% Courbes descripteur par descripteur, une quantite par panneau, legende
% partagee sous les panneaux.
n = size(G,1);
col = lines(K); mk = {'o','s','^','v','d','p'};
fig = figure;
t = tiledlayout(fig, 1, n, 'TileSpacing','compact', 'Padding','compact');
h = gobjects(K,1);
for q = 1:n
    D = G{q,1}; lim = G{q,3}; thr = G{q,4};
    ax = nexttile(t); hold(ax,'on'); box(ax,'on');
    Y = transformer(D, lim, St);
    for k = 1:K
        p = plot(ax, 1:S, Y(:,k), ['-' mk{k}], 'Color',col(k,:), 'LineWidth',1.0, ...
            'MarkerSize',3.5, 'MarkerFaceColor',col(k,:));
        if q == 1, h(k) = p; end
    end
    if ~isempty(thr), yline(ax, thr, '--k', 'LineWidth',1.0); end
    set(ax,'XTick',1:2:S,'XTickLabelRotation',0); xlim(ax, [0.7 S+0.3]);
    xlabel(ax, 'Scale level');
    if ischar(lim)
        ylim(ax, [0 St.plafond + 0.5]);
    else
        ylim(ax, lim);
    end
    title(ax, G{q,2}, 'FontWeight','bold');
end
lg = legend(h, names, 'Orientation','horizontal', 'Box','off');
lg.Layout.Tile = 'south';
paper_figstyle(fig, 15, 8);
exportgraphics(fig, fname, 'Resolution', St.dpi); close(fig);
end
