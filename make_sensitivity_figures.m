% version v10, 20/09/2026 : legendes remplacees par des intitules d'axes colores,
%                           panneaux resserres
% version v9, 20/09/2026 : trois panneaux superposes sur toute la largeur,
%                          etiquettes de configuration sous le dernier seulement,
%                          conventions reportees dans la legende de la figure
% version v8, 19/09/2026 : tri sans ambiguite, axe droit porte a 0.30, ordonnee
%                          du premier panneau nommee, repere des paires arretees
%                          remonte a mi-hauteur
% version v7, 19/09/2026 : axe propre pour la precision relative, conventions
%                          sur deux lignes
% version v6, 19/09/2026 : legendes reportees dans le quatrieme panneau, groupees
%                          par panneau ; precision relative a la place de C
% version v5, 19/09/2026 : une legende par panneau, le quatrieme ne portant
%                          plus que les conventions
% version v4, 19/09/2026 : barres conservees, disposition 2 x 2, canevas a la
%                          taille d'impression
% version v1 : barres en 1 x 3 sur un canevas de 34 cm
function make_sensitivity_figures(csvfile)
%MAKE_SENSITIVITY_FIGURES  Sensibilite des resultats au plan d'experience.
%   make_sensitivity_figures                     % defaut : volet2_results_paper.csv
%   make_sensitivity_figures('mon_fichier.csv')
%
% Produit : Fig_sensitivity_analysis.png (MEDIUM), Fig_sensitivity_low.png,
%           Fig_sensitivity_high.png
%
% CE QUI CHANGE PAR RAPPORT A LA v1. Le trace en barres par configuration est
% conserve. Deux causes d'illisibilite sont corrigees. Le canevas etait cree a
% 34 cm de large puis insere sur 15 cm, ce qui divisait toutes les tailles de
% police par 2,3 : il est desormais dessine a sa largeur d'impression. Et les
% trois panneaux etaient poses en 1 x 3, soit 4,3 cm pour 27 barres : ils sont
% en 2 x 2, ce qui porte chaque panneau a environ 7 cm, le quatrieme accueillant
% la legende et les conventions communes.

if nargin<1, csvfile='volet2_results_paper.csv'; end
T = readtable(csvfile);

% ===================== STYLE (tout se regle ici) =============================
S.dpi      = '-r600';
S.col1     = [0.23 0.46 0.69];   % bleu   (Power / R2 / FP)
S.col2     = [0.94 0.53 0.21];   % orange (ICC / C)
S.barw     = 0.42;
S.ylim12   = [0 1.05];
S.ylim3    = [0 0.15];
S.ylim_rp  = [0 0.30];           % axe droit ; la precision atteint 0.268 sous bruit fort
S.w_cm     = 15;                 % largeur d'impression
S.h_cm     = 13;                 % trois rangees resserrees
S.fs_tick  = 8;                  % etiquettes de configuration, verticales
S.fs_panel = 10;
S.fs_leg   = 9;
S.colx     = [0.55 0.12 0.12];   % repere des configurations arretees en amont
% =============================================================================

make_one(T, S, 'medium','Fig_sensitivity_analysis.png');
make_one(T, S, 'low',   'Fig_sensitivity_low.png');
make_one(T, S, 'high',  'Fig_sensitivity_high.png');
fprintf('Figures sauvees.\n');
end

% -----------------------------------------------------------------------------
function make_one(T, S, noise, fname)
% Trois panneaux superposes, chacun sur toute la largeur. Les configurations
% sont dans le meme ordre dans les trois : leurs etiquettes ne sont ecrites
% que sous le dernier. Les deux couleurs designent la premiere et la seconde
% barre de chaque paire ; elles sont reprises dans les intitules d'ordonnee,
% qui tiennent lieu de legende.
sel = T(strcmp(string(T.noise), noise), :);
% Tri par effort total, puis par G et par n a effort egal.
[~, ord] = sortrows([sel.G .* sel.n .* sel.R, sel.G, sel.n]);
sel = sel(ord, :);
x   = 1:height(sel);
labels = compose('G%d n%d R%d', sel.G, sel.n, sel.R);

f = figure;
t = tiledlayout(f, 3, 1, 'TileSpacing','tight', 'Padding','compact');

% ---- (1) Discrimination et fiabilite ----
ax = nexttile(t); hold(ax,'on');
bar(ax, x-S.barw/2, sel.power_F,  S.barw, 'FaceColor',S.col1);
bar(ax, x+S.barw/2, sel.mean_ICC, S.barw, 'FaceColor',S.col2);
yline(ax, 0.8, '--', 'Color',S.col1, 'HandleVisibility','off');
yline(ax, 0.6, '--', 'Color',S.col2, 'HandleVisibility','off');
% Pas de legende : l'intitule d'ordonnee porte les deux couleurs, comme les
% deux axes du panneau de correlation.
ylabel(ax, sprintf('{\\color[rgb]{%g,%g,%g}Power}, {\\color[rgb]{%g,%g,%g}ICC}', ...
       S.col1, S.col2), 'Interpreter','tex');
ylim(ax,S.ylim12); decor(ax, x, {}, S);

% ---- (2) Correlation, deux axes ----
ax = nexttile(t); hold(ax,'on');
yyaxis(ax,'left');
bar(ax, x-S.barw/2, sel.R2_recovery, S.barw, 'FaceColor',S.col1, 'EdgeColor','k');
yline(ax, 0.3, '--', 'Color',S.col1, 'HandleVisibility','off');
ylabel(ax,'Marginal R^2'); ylim(ax,S.ylim12); set(ax,'YColor',S.col1);
yyaxis(ax,'right');
bar(ax, x+S.barw/2, sel.rel_precision, S.barw, 'FaceColor',S.col2, 'EdgeColor','k');
ylabel(ax,'Relative precision'); ylim(ax,S.ylim_rp); set(ax,'YColor',S.col2);
yyaxis(ax,'left');
miss = isnan(sel.R2_recovery);           % paire arretee avant l'etage
if any(miss)
    plot(ax, x(miss), 0.5*ones(1,nnz(miss)), 'x', 'Color',S.colx, ...
         'MarkerSize',8, 'LineWidth',1.4, 'HandleVisibility','off');
end
decor(ax, x, {}, S);

% ---- (3) Erreur de type I ----
ax = nexttile(t); hold(ax,'on');
bar(ax, x, sel.false_positive_rate, 0.6, 'FaceColor',S.col1);
yline(ax, 0.05, '--r', 'HandleVisibility','off');
ylabel(ax,'False positive rate'); ylim(ax,S.ylim3);
decor(ax, x, labels, S);

paper_figstyle(f, S.w_cm, S.h_cm);
for a = findall(f,'Type','axes').'
    set(a, 'FontSize', S.fs_tick);       % apres le style commun
end
print(f, fname, '-dpng', S.dpi);
close(f);
end

% -----------------------------------------------------------------------------
function decor(ax, x, labels, S)
if isempty(labels)
    set(ax, 'XTick',x, 'XTickLabel',[]);
else
    set(ax, 'XTick',x, 'XTickLabel',labels, 'XTickLabelRotation',90, ...
            'TickLabelInterpreter','none');
end
xlim(ax, [0.4 numel(x)+0.6]);
box(ax,'on'); grid(ax,'on');
end
