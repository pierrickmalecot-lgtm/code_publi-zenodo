% version v2, 27/09/2026 : correction n/(n-1) de la variance bootstrap (NN1),
%                         comme dans paper_pipeline
% version v1, 02/09/2026 : stabilite des bornes bootstrap en fonction de B
%RUN_BCA_STABILITY  L'intervalle BCa tient-il a B = 1000.
%
%   PROBLEME. La correction de biais de BCa deplace les quantiles demandes vers
%   la queue de la distribution bootstrap. Pour une correction z0 de 0.84,
%   observee sur les descripteurs qui ne portent pas la relation, le quantile
%   haut demande vaut 0.99986 : sur B = 1000 tirages, la borne haute est le
%   maximum de l'echantillon, donc une seule realisation extreme. Une borne
%   definie par un seul tirage n'a aucune stabilite.
%
%   PRINCIPE DE MESURE. Le bootstrap est lance UNE FOIS avec NBLOCK * BLK
%   tirages par paire, puis decoupe en blocs independants. Chaque bloc, ou
%   groupe de blocs, fournit une execution independante a un B donne : 60 blocs
%   de 1000, 12 groupes de 5000, 3 de 20000, 1 de 60000. La dispersion des
%   bornes entre executions de meme B mesure directement leur instabilite,
%   sans qu'aucun ajustement soit refait.
%
%   TROIS INTERVALLES COMPARES, sur les memes tirages.
%     percentile   quantiles bruts alpha/2 et 1-alpha/2
%     de base      percentile inverse, 2*b1 - q_{1-alpha/2} et 2*b1 - q_{alpha/2}
%                  corrige un biais de translation sans demander de quantile
%                  extreme
%     BCa          corrige biais et asymetrie, au prix de quantiles extremes
%
%   LECTURE. Pour chaque paire et chaque B, l'etendue des bornes sur les
%   executions independantes, rapportee a l'ecart-type bootstrap. Une borne
%   utilisable doit avoir une etendue faible devant cet ecart-type. Si celle de
%   BCa reste grande la ou celle de l'intervalle de base est faible, BCa est a
%   ecarter a ce budget de calcul.
%
%   COUT. NPAIR * NBLOCK * BLK ajustements, parallelises par bloc. Avec les
%   valeurs par defaut, 240 000 ajustements.
%
%   Necessite : gen_volet1, paper_pipeline, paper_params.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

BLK    = 1000;    % tirages par bloc
NBLOCK = 60;      % nombre de blocs
ALPHA  = 0.05;
BGRID  = [1000 5000 20000 60000];   % doit diviser BLK*NBLOCK
NN1    = true;    % variance bootstrap multipliee par n/(n-1), comme le pipeline

P = paper_params();
rs = RandStream('Threefry','Seed',40007);      % tirage de reference, run_volet1
[X,Y] = gen_volet1(P,'nominal',rs);

% Quatre paires couvrant la plage de correction de biais observee : deux
% fortement biaisees, une intermediaire, et la paire porteuse s* = 6 sur D1,
% a partir de laquelle la reponse est construite.
PAIRS = [5 1; 6 2; 7 1; 6 1];
LBL   = {'D1 s=5','D2 s=6','D1 s=7','D1 s=6 (porteuse)'};
NPAIR = size(PAIRS,1);

% --- estimations ponctuelles et jackknife, une fois par paire -----------------
b1full = zeros(NPAIR,1); accel = zeros(NPAIR,1); separ = zeros(NPAIR,1);
for j = 1:NPAIR
    s = PAIRS(j,1); k = PAIRS(j,2);
    Xsk = X(:,:,s,k);
    [b1full(j), separ(j)] = fit_b1(Xsk, Y, mean(Xsk(:)), std(Xsk(:)));
    accel(j) = jack_accel(Xsk, Y);
    fprintf('%-20s  beta1 = %6.3f   se_par = %5.3f   acceleration = %+6.4f\n', ...
        LBL{j}, b1full(j), separ(j), accel(j));
end

% --- bootstrap par blocs ------------------------------------------------------
fprintf('\n=== Bootstrap : %d paires x %d blocs x %d tirages ===\n', ...
    NPAIR, NBLOCK, BLK);
tt = tic;
task = NPAIR*NBLOCK;
BS = zeros(BLK, task);
parfor t = 1:task
    j = 1 + floor((t-1)/NBLOCK);
    b = 1 + mod(t-1, NBLOCK);
    s = PAIRS(j,1); k = PAIRS(j,2);
    Xsk = X(:,:,s,k);
    mu_x = mean(Xsk(:)); sd_x = std(Xsk(:));
    rsb = RandStream('Threefry','Seed',200000 + 1000*j + b);
    BS(:,t) = boot_block(Xsk, Y, mu_x, sd_x, BLK, P, rsb);
end
fprintf('Termine en %.1f min.\n', toc(tt)/60);

% --- correction n/(n-1) -------------------------------------------------------
% Chaque pente est ecartee de l'estimation ponctuelle de sqrt(n/(n-1)), comme
% dans paper_pipeline. z0 est inchange ; les etendues rapportees a se_boot
% aussi, numerateur et denominateur etant multiplies par le meme facteur.
if NN1, fk = sqrt(P.n/(P.n-1)); else, fk = 1; end
for j = 1:NPAIR
    cols = (j-1)*NBLOCK + (1:NBLOCK);
    BS(:,cols) = b1full(j) + fk*(BS(:,cols) - b1full(j));
end

% --- analyse ------------------------------------------------------------------
for j = 1:NPAIR
    cols = (j-1)*NBLOCK + (1:NBLOCK);
    v_all = BS(:,cols); v_all = v_all(:);
    v_all = v_all(~isnan(v_all));
    seb = std(v_all);
    z0_all = norminv(clampp(mean(v_all < b1full(j)), numel(v_all)));

    fprintf('\n================ %s ================\n', LBL{j});
    fprintf('beta1 = %.3f | se_boot = %.3f | se_boot/se_par = %.2f | ', ...
        b1full(j), seb, seb/separ(j));
    fprintf('z0 = %.3f | accel = %+.4f\n', z0_all, accel(j));
    a2 = bca_alpha(z0_all, accel(j), norminv(1-ALPHA/2));
    fprintf('Quantile haut demande par BCa : %.5f\n', a2);
    for B = BGRID
        if ~isnan(a2)
            fprintf('  B = %6d  ->  rang %.1f sur %d, %.1f tirages au-dela\n', ...
                B, a2*B, B, (1-a2)*B);
        end
    end

    fprintf('\n  methode      B      borne basse moyenne   etendue/se   borne haute moyenne   etendue/se\n');
    for B = BGRID
        nrun = BLK*NBLOCK/B;
        if nrun < 2, continue; end
        per = B/BLK;                      % blocs par execution
        lo = zeros(nrun,3); hi = zeros(nrun,3);
        for r = 1:nrun
            cc = cols((r-1)*per + (1:per));
            v = BS(:,cc); v = v(:); v = v(~isnan(v));
            [lo(r,:), hi(r,:)] = three_ci(v, b1full(j), accel(j), ALPHA);
        end
        nm = {'percentile','de base   ','BCa       '};
        for meth = 1:3
            fprintf('  %s %6d  %19.3f  %11.2f  %20.3f  %11.2f\n', ...
                nm{meth}, B, mean(lo(:,meth)), range(lo(:,meth))/seb, ...
                mean(hi(:,meth)), range(hi(:,meth))/seb);
        end
    end
end

save('ws_run_bca_stability.mat','BS','PAIRS','LBL','b1full','accel','separ', ...
     'BLK','NBLOCK','BGRID','ALPHA','P','NN1');
fprintf('\nEspace de travail enregistre dans ws_run_bca_stability.mat\n');

% =============================================================================
function v = boot_block(Xsk, Y, mu_x, sd_x, nb, P, rsb)
% Un bloc de nb tirages bootstrap specimen, stratifie par condition.
[N,~] = size(Xsk); G = P.G; n = P.n;
v = nan(nb,1);
for b = 1:nb
    idx = zeros(N,1);
    for g = 1:G
        base = (g-1)*n;
        idx(base+(1:n)) = base + randi(rsb, n, n, 1);
    end
    try
        v(b) = fit_b1(Xsk(idx,:), Y(idx,:), mu_x, sd_x);
    catch
    end
end
end
% -----------------------------------------------------------------------------
function [b1, se1] = fit_b1(Xsk, Y, mu_x, sd_x)
% Meme ajustement que lme_b1_only de paper_pipeline : REML, standardisation
% fixee sur l'echantillon complet.
[N,R] = size(Xsk);
x = (Xsk(:)-mu_x)/sd_x; y = Y(:);
T = table(y, x, nominal(repelem((1:N)',R,1)), ...
          'VariableNames',{'Y','X','Sample'});
lme = fitlme(T,'Y ~ X + (1|Sample)','FitMethod','REML');
b1  = lme.Coefficients.Estimate(2);
se1 = lme.Coefficients.SE(2);
end
% -----------------------------------------------------------------------------
function a = jack_accel(Xsk, Y)
% Acceleration BCa par jackknife au niveau specimen.
[N,~] = size(Xsk);
mu_x = mean(Xsk(:)); sd_x = std(Xsk(:));
th = nan(N,1);
for i = 1:N
    keep = true(N,1); keep(i) = false;
    try, th(i) = fit_b1(Xsk(keep,:), Y(keep,:), mu_x, sd_x); catch, end
end
th = th(~isnan(th));
if numel(th) < 3, a = 0; return; end
d = mean(th) - th; den = 6*(sum(d.^2))^1.5;
if den < eps, a = 0; else, a = sum(d.^3)/den; end
end
% -----------------------------------------------------------------------------
function [lo, hi] = three_ci(v, b1, accel, alpha)
% Renvoie les bornes des trois intervalles : percentile, de base, BCa.
q = quantile(v, [alpha/2, 1-alpha/2]);
lo = zeros(1,3); hi = zeros(1,3);
lo(1) = q(1);            hi(1) = q(2);              % percentile
lo(2) = 2*b1 - q(2);     hi(2) = 2*b1 - q(1);       % de base
z0 = norminv(clampp(mean(v < b1), numel(v)));
a1 = bca_alpha(z0, accel, norminv(alpha/2));
a2 = bca_alpha(z0, accel, norminv(1-alpha/2));
if isnan(a1) || isnan(a2) || a1 >= a2
    lo(3) = q(1); hi(3) = q(2);
else
    qb = quantile(v, [a1 a2]);
    lo(3) = qb(1); hi(3) = qb(2);
end
end
% -----------------------------------------------------------------------------
function p = clampp(p, nb)
% Borne la proportion pour que z0 reste fini.
p = min(max(p, 1/(2*nb)), 1 - 1/(2*nb));
end
% -----------------------------------------------------------------------------
function a = bca_alpha(z0, accel, z)
den = 1 - accel*(z0 + z);
if den <= 0, a = NaN; else, a = normcdf(z0 + (z0 + z)/den); end
end
