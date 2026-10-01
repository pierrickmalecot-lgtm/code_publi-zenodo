% version v9, 19/09/2026 : une seule barre de couleur pour les deux matrices BC
% version v1, 05/09/2026 : regeneration des figures sans recalcul
%MAKE_ALL_FIGURES  Retrace toutes les figures du papier a partir des espaces de
%   travail deja enregistres, au format d'impression.
%
%   POURQUOI. Jusqu'ici chaque figure etait tracee a la fin du script de calcul
%   qui l'alimente. Changer une mise en forme imposait donc de relancer le
%   calcul : vingt-six minutes pour run_boot_vs_param, une heure pour le
%   dimensionnement. Or les donnees sont deja dans les fichiers ws_*.mat.
%   Ce script les relit et redessine, en quelques secondes.
%
%   FORMAT. Toutes les figures sont composees pour 15 cm de large, la valeur
%   d'une figure sur deux colonnes chez l'editeur, avec des polices de 9 pt.
%   La mise en forme est appliquee par paper_figstyle, seul endroit a modifier
%   si le format change.
%
%   IMPORTANT. Les figures etant desormais dessinees a leur taille finale, elles
%   doivent etre incluses dans le document SANS mise a l'echelle, par
%   \includegraphics[width=15cm]{...} et non width=\textwidth, faute de quoi
%   LaTeX annule le travail fait ici.
%
%   PREREQUIS. Les fichiers ws_*.mat correspondants doivent exister, donc les
%   etudes doivent avoir tourne au moins une fois. Une etude absente est
%   signalee et sautee, les autres sont produites.
%
%   USAGE.
%       make_all_figures          % toutes celles dont les donnees existent
%       make_all_figures('perm_size')   % une seule
function make_all_figures(varargin)

figs = { ...
  'boot_vs_param',   'ws_run_boot_vs_param.mat',   @fig_boot_vs_param   ; ...
  'perm_size',       'ws_run_perm_size.mat',       @fig_perm_size       ; ...
  'curvature_noise', 'ws_run_curvature_noise.mat', @fig_curvature_noise ; ...
  'loso_sweep',      'ws_run_loso_sweep.mat',      @fig_loso_sweep      ; ...
  'randslope_R',     'ws_run_randslope_R.mat',     @fig_randslope_R     ; ...
  'decision_matrix', 'ws_run_volet1.mat',          @fig_decision_matrix ; ...
  'q_indicator',     'ws_run_volet1.mat',          @fig_q_indicator     ; ...
  'contrast_sweep',  'ws_run_contrast_sweep.mat',  @fig_contrast_sweep  ; ...
  'scenarios_bc',    'ws_run_scenarios.mat',       @fig_scenarios_bc    ; ...
  'unreliable',      'ws_make_unreliable_figure.mat', @fig_unreliable   ; ...
};

if nargin >= 1
    keep = strcmpi(figs(:,1), varargin{1});
    if ~any(keep)
        error('make_all_figures:unknown', ...
              'Figure inconnue : %s. Disponibles : %s.', varargin{1}, ...
              strjoin(figs(:,1)', ', '));
    end
    figs = figs(keep,:);
end

fprintf('=== Regeneration des figures, format 15 cm ===\n\n');
nOK = 0; nSkip = 0;
for i = 1:size(figs,1)
    nom = figs{i,1}; wsf = figs{i,2}; fun = figs{i,3};
    if exist(wsf,'file') ~= 2
        fprintf('  saute   %-18s (%s absent)\n', nom, wsf);
        nSkip = nSkip + 1;
        continue
    end
    D = load(wsf);
    try
        out = fun(D);
        fprintf('  ok      %-18s -> %s\n', nom, out);
        nOK = nOK + 1;
    catch ME
        fprintf('  ECHEC   %-18s : %s\n', nom, ME.message);
    end
end
fprintf('\n%d figure(s) produite(s), %d sautee(s).\n', nOK, nSkip);
fprintf(['Inclure dans le document par \\includegraphics[width=15cm]{...},\n' ...
         'sans mise a l echelle.\n']);
end

% =============================================================================
function nom = fig_boot_vs_param(D)
% Rapport bootstrap/parametrique et biais normalise contre echelle, D1 a D3.
K = D.P.K;
kShow = [];
for k = 1:K
    if any(D.freq(:,k) >= D.FREQ_MIN), kShow(end+1) = k; end %#ok<AGROW>
end
col = lines(numel(kShow));
f = figure;
for panel = 1:2
    subplot(1,2,panel); hold on;
    for j = 1:numel(kShow)
        k = kShow(j);
        ss = find(D.freq(:,k) >= D.FREQ_MIN)';
        med = nan(size(ss)); lo = med; hi = med;
        for a = 1:numel(ss)
            if panel == 1, v = D.ratio(:,ss(a),k); else, v = D.biasN(:,ss(a),k); end
            v = v(~isnan(v));
            med(a) = median(v); q = quantile(v,[0.25 0.75]);
            lo(a) = q(1); hi(a) = q(2);
        end
        fill([ss fliplr(ss)], [lo fliplr(hi)], col(j,:), ...
             'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
        plot(ss, med, '-o', 'Color',col(j,:), 'LineWidth',1.2, 'MarkerSize',4);
    end
    xlabel('Scale index s');
    if panel == 1
        ylabel('se_{bootstrap} / se_{parametric}');
        yline(1,':','HandleVisibility','off');
        legend(D.names(kShow),'Location','northeast');
    else
        ylabel('Normalised bootstrap bias');
        yline(0,':','HandleVisibility','off');
    end
    grid on; box on;
end
paper_figstyle(f, 15, 7);
nom = 'Fig_boot_vs_param.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_perm_size(D)
% Taux d'erreur de type I contre nombre de specimens, D4 et D5.
LEG = {'branching','permutation (F)','permutation (Welch)'};
mk  = {'-o','-s','-^'};
tit = {'D4: pure region noise','D5: specimen + region noise'};
f = figure;
for panel = 1:2
    subplot(1,2,panel); hold on;
    for im = 1:numel(LEG)
        plot(D.NGRID, D.alpha_null(:,im,panel), mk{im}, ...
             'LineWidth',1.2, 'MarkerSize',5);
    end
    yline(0.05,'--','HandleVisibility','off');
    set(gca,'XScale','log','XTick',D.NGRID,'XTickLabel',string(D.NGRID));
    xlabel('Specimens per condition, n');
    ylabel('Type-I error rate');
    title(tit{panel});
    % La legende est placee hors des courbes et le haut de l'axe degage : a
    % 10 pt elle mordait sur la courbe superieure du panneau gauche.
    if panel == 1
        yl = ylim; ylim([yl(1) yl(1) + 1.30*(yl(2)-yl(1))]);
        legend(LEG,'Location','northeast');
    end
    grid on; box on;
end
paper_figstyle(f, 15, 7.5);
nom = 'Fig_perm_size.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_curvature_noise(D)
% Degradation du R2 marginal et du Spearman avec la courbure, par niveau de bruit.
mk = {'-o','-s','-^'};
nN = numel(D.noiseLbl);
f = figure;
subplot(1,2,1); hold on;
for j = 1:nN
    plot(D.curvGrid, D.R2(:,j), mk{j}, 'LineWidth',1.2,'MarkerSize',4);
end
xlabel('Curvature rate c'); ylabel('R^2_{marginal}');
legend(D.noiseLbl,'Location','southwest');
grid on; box on; ylim([0 1]); title('Linear fit');
subplot(1,2,2); hold on;
for j = 1:nN
    plot(D.curvGrid, D.RHO(:,j), mk{j}, 'LineWidth',1.2,'MarkerSize',4);
end
xlabel('Curvature rate c'); ylabel('Spearman \rho');
grid on; box on; ylim([0 1]); title('Rank association');
paper_figstyle(f, 15, 7);
nom = 'Fig_curvature_noise.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_loso_sweep(D)
% Frequence du drapeau fragile contre amplitude de l'anomalie.
f = figure;
plot(D.ampGrid, D.freqTgt, '-o', 'LineWidth',1.2, 'MarkerSize',4); hold on;
plot(D.ampGrid, D.freqAny, '--s', 'LineWidth',1.0, 'MarkerSize',4);
yline(0.5,':','HandleVisibility','off');
xlabel('ROI dispersion inflation factor of the anomalous specimen');
ylabel('Frequency of the fragility flag');
legend({'target pair','any pair'},'Location','southeast');
ylim([0 1.05]); xlim([min(D.ampGrid) max(D.ampGrid)]);
grid on; box on;
paper_figstyle(f, 15, 7.5);
nom = 'Fig_loso_sweep.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_randslope_R(D)
% Selection du modele a pente aleatoire : balayage en R contre balayage en n.
f = figure;
subplot(1,2,1); hold on;
plot(D.GRID, D.res(1).pA, '-o', 'LineWidth',1.2, 'MarkerSize',4);
plot(D.GRID, D.res(2).pA, '--s','LineWidth',1.2, 'MarkerSize',4);
yline(0.5,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',D.GRID,'XTickLabel',string(D.GRID));
xlabel('Regions per specimen, or specimens per condition');
ylabel('P(AIC selects random slope)');
legend({'sweep in R','sweep in n'},'Location','northwest');
grid on; box on; ylim([0 1]);
subplot(1,2,2); hold on;
plot(D.GRID, D.res(1).dA, '-o', 'LineWidth',1.2, 'MarkerSize',4);
plot(D.GRID, D.res(2).dA, '--s','LineWidth',1.2, 'MarkerSize',4);
yline(0,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',D.GRID,'XTickLabel',string(D.GRID));
xlabel('Regions per specimen, or specimens per condition');
ylabel('\DeltaAIC (reduced - full)');
grid on; box on;
paper_figstyle(f, 15, 7);
nom = 'Fig_randslope_R.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_decision_matrix(D)
% Matrice de decision des trois portes.
%
% CORRECTION. La version precedente posait 'YTick',1:6 avec six etiquettes de
% descripteurs alors que le panel n'en compte que cinq, la sixieme entree de
% descriptorNames datant d'un panel anterieur. L'etiquetage etait donc decale
% d'une ligne. Le nombre de lignes est desormais lu sur les donnees.
K   = D.P.K;
dn  = cellstr(string(D.descriptorNames(1:K)));
f = figure;
imagesc(D.out.decision');
colormap(gca,[0.90 0.90 0.90; 0.70 0.85 1.00; 0.30 0.65 0.90; 0.15 0.70 0.30]);
clim([-0.25 3.25]);
colorbar('Ticks',[0.375 1.125 1.875 2.625],'TickLabels', ...
    {'Excluded','Discrim. only','Discrim. + reliable','Retained'});
set(gca,'YTick',1:K,'YTickLabel',dn,'XTick',1:D.P.S);
xlabel('Scale level'); ylabel('Descriptor');
title('Three-stage filtering decision matrix');
paper_figstyle(f, 15, 7);
nom = 'Fig_decision_matrix.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_q_indicator(D)
% Indicateur composite Q par paire descripteur-echelle.
%
% CORRECTION. Q etait alloue en P.S par 6 alors que la boucle n'en remplit que
% cinq colonnes, le panel comptant cinq descripteurs. La sixieme ligne de la
% carte etait donc uniformement nulle et s'affichait comme un descripteur
% supplementaire sans aucune valeur. Q est desormais alloue au nombre reel de
% descripteurs.
K = D.P.K; S = D.P.S;
dn = cellstr(string(D.descriptorNames(1:K)));
Q = zeros(S,K);
for s = 1:S
    for k = 1:K
        qd = min(D.out.eta2(s,k)/0.14, 1);
        qi = min(D.out.ICC(s,k)/D.P.icc_thr, 1);
        if isnan(D.out.C(s,k)), qb = 0; else, qb = min(D.out.C(s,k)/D.P.C_thr, 1); end
        Q(s,k) = (qd + qi + qb)/3;
    end
end
f = figure;
imagesc(Q'); colorbar; clim([0 1]); colormap(gca, parula);
set(gca,'YTick',1:K,'YTickLabel',dn,'XTick',1:S);
xlabel('Scale level'); ylabel('Descriptor');
title('Composite reliability indicator Q_{s,k}');
paper_figstyle(f, 15, 7);
nom = 'Fig_Q_indicator.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_contrast_sweep(D)
% Fiabilite contre contraste, et contre nombre de conditions a variance egalisee.
f = figure;
subplot(1,2,1); hold on;
plot(D.ratioEff, D.iccA, '-o','LineWidth',1.2,'MarkerSize',4);
cHat = mean( D.iccA(:)./(1-D.iccA(:)) ./ (D.ratioEff(:).^2) );
rr = linspace(min(D.ratioEff), max(D.ratioEff), 200);
plot(rr, cHat*(rr.^2)./(cHat*(rr.^2)+1), '--','LineWidth',1.0);
yline(D.P.icc_thr,':','HandleVisibility','off');
xlabel('Descriptor contrast  \Delta / \sigma_{within}');
ylabel('ICC'); ylim([0 1]); grid on; box on;
legend({'simulated','fitted'},'Location','southeast');
title('Reliability against contrast');
subplot(1,2,2); hold on;
plot(D.Ggrid, D.iccB, '-o','LineWidth',1.2,'MarkerSize',4);
plot(D.Ggrid, D.powB, '-s','LineWidth',1.2,'MarkerSize',4);
yline(D.P.icc_thr,':','HandleVisibility','off');
xlabel('Number of conditions G'); ylabel('Value'); ylim([0 1.05]);
xticks(D.Ggrid); grid on; box on;
legend({'ICC','Power'},'Location','southeast');
title('At matched between-condition variance');
paper_figstyle(f, 15, 7);
nom = 'Fig_contrast_sweep.png';
print(f, nom, '-dpng', '-r300');
end

% -----------------------------------------------------------------------------
function nom = fig_scenarios_bc(D)
% Matrices de decision des scenarios null et deux echelles.
%
% CORRECTION v1. La version precedente posait 'YTick',1:6 avec six etiquettes
% alors que le panel compte cinq descripteurs, la sixieme entree datant d'un
% panel anterieur. Le nombre de lignes est desormais lu sur les donnees.
%
% CORRECTION v2. Les deux panneaux partagent la meme palette et les memes
% bornes, et portaient chacun leur barre de couleur. Une seule barre sert les
% deux, ce qui rend aux panneaux la place qu'elles occupaient. Graduations
% impaires et cases proportionnees, comme dans les autres matrices.
K  = D.P.K;
dn = cellstr(string(D.dn(1:K)));
tt = cellstr(string(D.tt));
dd = {D.outB.decision, D.outC.decision};
f = figure;
t = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');
for pnl = 1:2
    ax = nexttile(t);
    imagesc(ax, dd{pnl}');
    colormap(ax, D.cmap); clim(ax, [-0.25 3.25]);
    set(ax,'YTick',1:K,'YTickLabel',dn, ...
           'XTick',1:2:D.P.S,'XTickLabelRotation',0);
    pbaspect(ax, [D.P.S 2*K 1]);
    xlabel(ax, 'Scale level');
    if pnl == 1, ylabel(ax, 'Descriptor'); end
    title(ax, tt{pnl}, 'FontWeight','bold');
end
cb = colorbar('Ticks',[0.375 1.125 1.875 2.625],'TickLabels', ...
    {'Excluded','Discrim. only','Discrim.+reliable','Retained'});
cb.Layout.Tile = 'east';
paper_figstyle(f, 15, 8);
nom = 'Fig_decision_matrices_BC.png';
print(f, nom, '-dpng', '-r300');
end
% -----------------------------------------------------------------------------
function nom = fig_unreliable(D)
% Descripteur correle mais peu fiable : ICC, R2 marginal et Spearman contre
% echelle, compares au descripteur porteur.
mu = @(A) mean(A,1,'omitnan');
sc = 1:D.P.S;
f = figure;
subplot(1,3,1);
plot(sc, mu(D.ICCr), '-o','LineWidth',1.2,'MarkerSize',3.5); hold on;
plot(sc, mu(D.ICCu), '-s','LineWidth',1.2,'MarkerSize',3.5);
yline(D.P.icc_thr,'--','HandleVisibility','off');
xlabel('Scale level'); ylabel('ICC'); ylim([0 1]); grid on; box on;
legend({'D1','unreliable carrier'},'Location','northwest');
title('Reliability');
subplot(1,3,2);
plot(sc, mu(D.R2r), '-o','LineWidth',1.2,'MarkerSize',3.5); hold on;
plot(sc, mu(D.R2u), '-s','LineWidth',1.2,'MarkerSize',3.5);
yline(D.P.r2_thr,'--','HandleVisibility','off');
xlabel('Scale level'); ylabel('R^2_{marginal} (gates disabled)');
ylim([0 1]); grid on; box on; title('Correlation');
subplot(1,3,3);
plot(sc, mu(D.RHr), '-o','LineWidth',1.2,'MarkerSize',3.5); hold on;
plot(sc, mu(D.RHu), '-s','LineWidth',1.2,'MarkerSize',3.5);
xlabel('Scale level'); ylabel('Spearman \rho');
ylim([-0.2 1]); grid on; box on; title('Monotonicity');
paper_figstyle(f, 15, 6);
nom = 'Fig_unreliable_curves.png';
print(f, nom, '-dpng', '-r300');
end
