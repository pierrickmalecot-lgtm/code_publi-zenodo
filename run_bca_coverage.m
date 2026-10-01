% version v2, 27/09/2026 : le pipeline applique n/(n-1) par defaut ; le script
%                         compare BCa avec et sans correction a partir des
%                         champs b1_*_raw, sans calcul maison du facteur
% version v1, 17/09/2026 : couverture des intervalles de la pente, cas exact
%RUN_BCA_COVERAGE  L'intervalle BCa de la pente couvre-t-il la valeur vraie
%   avec la probabilite annoncee, a n = 5 specimens par condition.
%
%   CONTEXTE. En validant la correction d'attenuation, le bootstrap stratifie
%   par condition s'est montre trop peu disperse a n = 5 : les intervalles qui
%   en derivent couvraient la pente vraie dans 87 pour cent des cas environ,
%   pour 95 vises. Le bootstrap du pipeline (boot_C) est stratifie de la meme
%   facon. Un controle approche en Python, moindres carres groupes a la place
%   du modele mixte, donne 90 pour cent pour BCa et 91 pour le percentile. Ce
%   script mesure la couverture sur le pipeline lui-meme.
%
%   PROTOCOLE. Scenario nominal, paire D1 a s*, la seule dont la pente vraie
%   est connue : Y y est construit directement sur X(i,r,s*,1) avec la pente
%   P.beta1. Le pipeline standardise X sur l'echantillon complet ; la pente
%   vraie dans cette unite vaut donc P.beta1 * sd(X), recalculee sur chaque
%   replicat. Budget bootstrap du papier, B = 1000.
%
%   INTERVALLES COMPARES, a 95 pour cent.
%     BCa corrige      celui du papier : variance x n/(n-1), puis BCa
%     BCa brut         meme bootstrap, sans la correction n/(n-1)
%     percentile brut  quantiles du bootstrap sans correction
%     parametrique     b1 +/- z * se de fitlme, echantillon de specimens fixe
%     normal brut      b1 +/- z * ecart-type bootstrap sans correction
%     normal corrige   b1 +/- z * ecart-type bootstrap corrige
%   RESULTAT v1 (17/09/2026) : BCa brut 91.5 %, parametrique 75.0 %, normal
%   corrige 95.5 %. La couverture de BCa corrige n'avait pas ete mesuree.
%
%   LECTURE. Avec 200 replicats, l'erreur type sur une couverture de 95 pour
%   cent vaut 1.5 point. Un intervalle a 90 pour cent ou moins est donc trop
%   etroit de facon mesurable. Le nombre de replicats ou D1 a s* n'est pas
%   retenu au verrou 3 est rapporte ; ces replicats sont exclus.
%
%   GRAINES. Niveau 0 de paper_seed : memes jeux que les autres etudes
%   nominales sans balayage. Le flux bootstrap est celui du travailleur, comme
%   dans run_boot_vs_param.
%
%   COUT. NREP * (B + N) ajustements mixtes, environ 205 000 a NREP = 200.
%   Compter de 15 a 30 minutes selon le nombre de travailleurs.
%
%   Necessite : gen_volet1, paper_pipeline, paper_params, paper_seed.
clear; clc;
pool = gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP   = 200;
B_BOOT = 1000;
ALPHA  = 0.05;

P  = paper_params();
s0 = P.s_star; k0 = 1;
z  = norminv(1 - ALPHA/2);
NAMES = {'BCa corrige','BCa brut','percentile brut','parametrique', ...
         'normal brut','normal corrige'};
nI = numel(NAMES);

cov   = nan(NREP, nI);
wid   = nan(NREP, nI);
isBCa = nan(NREP, 1);
z0v   = nan(NREP, 1);

fprintf('=== Couverture de l intervalle de la pente, D1 a s = %d ===\n', s0);
fprintf('%d replicats, B = %d, n = %d specimens par condition\n\n', NREP, B_BOOT, P.n);
tt = tic;

parfor r = 1:NREP
    rs = RandStream('Threefry','Seed', paper_seed('nominal', r));
    [X, Y] = gen_volet1(P, 'nominal', rs);
    o = paper_pipeline(X, Y, P, struct('bootstrap', true, 'B', B_BOOT, ...
                                       'bootPairs', [s0 k0], 'nn1', true));
    if o.decision(s0,k0) < 3 || isnan(o.b1_lo(s0,k0)), continue; end

    xs  = X(:,:,s0,k0);
    bt  = P.beta1 * std(xs(:));             % pente vraie, unite standardisee
    b1  = o.beta1(s0,k0);
    sb  = o.b1_se_raw(s0,k0);               % ecart-type bootstrap sans correction
    sc  = o.b1_se(s0,k0);                   % ecart-type bootstrap corrige
    sp  = o.b1_se_par(s0,k0);
    lo  = [o.b1_lo(s0,k0), o.b1_lo_raw(s0,k0), o.b1_pct_lo_raw(s0,k0), ...
           b1 - z*sp, b1 - z*sb, b1 - z*sc];
    hi  = [o.b1_hi(s0,k0), o.b1_hi_raw(s0,k0), o.b1_pct_hi_raw(s0,k0), ...
           b1 + z*sp, b1 + z*sb, b1 + z*sc];

    cov(r,:)  = double(lo <= bt & bt <= hi);
    wid(r,:)  = (hi - lo) / abs(b1);
    isBCa(r)  = double(strcmp(o.b1_method{s0,k0}, 'BCa'));
    z0v(r)    = o.b1_z0(s0,k0);
end
fprintf('Termine en %.1f min.\n\n', toc(tt)/60);

v = ~isnan(cov(:,1));
fprintf('Replicats retenus : %d sur %d\n', nnz(v), NREP);
fprintf('BCa effectivement applique : %.1f %%, z0 median %.3f\n\n', ...
    100*mean(isBCa(v)), median(z0v(v)));

fprintf('%-16s %12s %12s\n', 'intervalle', 'couverture', 'larg. rel.');
for j = 1:nI
    fprintf('%-16s %11.1f %% %12.3f\n', NAMES{j}, 100*mean(cov(v,j)), median(wid(v,j)));
end

fprintf(['\nLECTURE. La couverture visee est 95 %%, a 1.5 point pres. La largeur\n' ...
         'relative est la largeur rapportee a la pente estimee. La comparaison\n' ...
         'utile est BCa corrige contre BCa brut : elle mesure ce que la\n' ...
         'correction n/(n-1) apporte a l intervalle publie.\n']);

save('ws_run_bca_coverage.mat','cov','wid','isBCa','z0v','NAMES','NREP','B_BOOT','P');
fprintf('\nws_run_bca_coverage.mat ecrit.\n');
