% version v9, 19/09/2026 : bandes interquartiles, bande injectee marquee,
%                          legende sur les trois panneaux, retrace sans recalcul
% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/08/2026 : figure du descripteur correle mais peu fiable (4.9)
%MAKE_UNRELIABLE_FIGURE  Trois panneaux, ICC / R2 marginal / Spearman rho
%   contre l'echelle, pour le descripteur correle mais peu fiable (DU) avec D1
%   en reference. Le point de la figure est que les deux courbes se superposent
%   sur la correlation et divergent sur la fiabilite : c'est ce que le texte de
%   la section 4.9 affirme sans le montrer.
%
%   Le R2 est calcule avec les verrous amont desactives, faute de quoi DU
%   n'atteindrait jamais l'etage de correlation.
%   RETRACE SANS RECALCUL. Les 200 replicats sont conserves en entier dans
%   ws_make_unreliable_figure.mat. Si ce fichier est present, la figure est
%   retracee a partir de lui ; le recalcul n'a lieu qu'en son absence, ou si
%   RECALC vaut true.
clear; clc;
RECALC = false;
WS     = 'ws_make_unreliable_figure.mat';

P = paper_params(); P.K = 6; kU = 6; kRef = 1; NREP = 200;

if ~RECALC && exist(WS,'file')
    load(WS,'ICCu','ICCr','R2u','R2r','RHu','RHr','NREP','P');
    fprintf('Donnees relues dans %s, %d replicats.\n', WS, NREP);
else
    pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
    Pn = P; Pn.q=1; Pn.eta2_thr=0; Pn.icc_thr=0; Pn.r2_thr=-Inf;   % verrous ouverts

    ICCu=nan(NREP,P.S); ICCr=nan(NREP,P.S);
    R2u =nan(NREP,P.S); R2r =nan(NREP,P.S);
    RHu =nan(NREP,P.S); RHr =nan(NREP,P.S);

    fprintf('=== Courbes DU vs D1, %d replicats ===\n', NREP);
    tt=tic;
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('unreliable', r));   % memes graines que run_unreliable
        [X,Y] = gen_volet1(P,'unreliable',rs);
        o  = paper_pipeline(X,Y,P);
        of = paper_pipeline(X,Y,Pn);
        ICCu(r,:)=o.ICC(:,kU)';  ICCr(r,:)=o.ICC(:,kRef)';
        R2u(r,:) =of.R2m(:,kU)'; R2r(r,:) =of.R2m(:,kRef)';
        RHu(r,:) =o.rho(:,kU)';  RHr(r,:) =o.rho(:,kRef)';
    end
    fprintf('Termine en %.1f min.\n', toc(tt)/60);
    save(WS,'ICCu','ICCr','R2u','R2r','RHu','RHr','NREP','P');
end

mu = @(A) mean(A,1,'omitnan');

f = figure;
t = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');

panneau(nexttile(t), ICCr, ICCu, 'ICC', [0 1], P.icc_thr, 'Reliability');
panneau(nexttile(t), R2r,  R2u,  'R^2_{marginal} (gates disabled)', [0 1], P.r2_thr, 'Correlation');
panneau(nexttile(t), RHr,  RHu,  'Spearman \rho', [-0.2 1], [], 'Monotonicity');

paper_figstyle(f, 15, 7);
print(f,'Fig_unreliable_curves.png','-dpng','-r300');
fprintf('Fig_unreliable_curves.png sauvee.\n');

% --- valeurs citables dans le texte ------------------------------------------
band = 5:7;
fprintf('\n=== Sur la bande 5-7 ===\n');
fprintf('  ICC   : D1 = %.3f | DU = %.3f\n', mean(mu(ICCr(:,band))), mean(mu(ICCu(:,band))));
fprintf('  R2    : D1 = %.3f | DU = %.3f\n', mean(mu(R2r(:,band))),  mean(mu(R2u(:,band))));
fprintf('  rho   : D1 = %.3f | DU = %.3f\n', mean(mu(RHr(:,band))),  mean(mu(RHu(:,band))));
fprintf('\nrho de DU hors bande (toutes echelles) : %.3f\n', mean(mu(RHu)));

% -----------------------------------------------------------------------------
function panneau(ax, Aref, Au, ylab, yl, thr, titre)
% Moyenne sur les replicats, entouree de l'intervalle interquartile. La bande
% ou la reponse fonctionnelle est construite est marquee en fond, de sorte que
% le lecteur voie ce que les courbes font a l'interieur et a l'exterieur.
mu  = @(A) mean(A,1,'omitnan');
q25 = @(A) quantile(A,0.25,1);
q75 = @(A) quantile(A,0.75,1);
S   = size(Aref,2); sc = 1:S; COL = lines(2); BAND = [5 7];
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
patch(ax, [BAND(1)-0.5 BAND(2)+0.5 BAND(2)+0.5 BAND(1)-0.5], ...
          [yl(1) yl(1) yl(2) yl(2)], [0.92 0.92 0.92], ...
          'EdgeColor','none','HandleVisibility','off');
A = {Aref, Au}; mkr = {'-o','-s'}; h = gobjects(2,1);
for i = 1:2
    fill(ax, [sc fliplr(sc)], [q25(A{i}) fliplr(q75(A{i}))], COL(i,:), ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    h(i) = plot(ax, sc, mu(A{i}), mkr{i}, 'Color',COL(i,:), ...
        'LineWidth',1.4,'MarkerSize',4,'MarkerFaceColor',COL(i,:));
end
if ~isempty(thr) && isfinite(thr)
    yline(ax, thr, '--', 'HandleVisibility','off');
end
set(ax,'XTick',1:2:S,'XTickLabelRotation',0);
xlabel(ax,'Scale level'); ylabel(ax,ylab); ylim(ax,yl); xlim(ax,[0.5 S+0.5]);
title(ax, titre, 'FontWeight','bold');
legend(ax, h, {'D1','unreliable carrier'}, 'Location','best', 'Box','off');
end
