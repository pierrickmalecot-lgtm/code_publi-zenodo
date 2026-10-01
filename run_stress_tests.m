% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1 : tests de robustesse aux violations d'hypotheses (section 6.4)
%RUN_STRESS_TESTS  Volet 1c : degradation du pipeline sous violation des
%   hypotheses statistiques, une hypothese a la fois, plan de reference.
%
%   Trois familles de violation :
%     (a) independance des ROI  : bruit intra AR(1) le long de r, rho_spat
%     (b) normalite             : bruit intra log-normal (scenario existant)
%     (c) linearite             : reponse saturante, taux de courbure P.curv
%
%   Necessite gen_volet1.m avec les scenarios 'autocorr' et 'curvature',
%   et paper_params.m avec P.rho_spat et P.curv.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P=paper_params(); NREP=200;

kRef=1; kNull=4; s_star=P.s_star; bandRef=5:7;

% --- definition des lignes du tableau ----------------------------------------
% {libelle, scenario, champ modifie, valeur}
rows = { ...
    'reference (nominal)'          , 'nominal'   , ''         , NaN  ; ...
    'AR(1) rho = 0.3'              , 'autocorr'  , 'rho_spat' , 0.3  ; ...
    'AR(1) rho = 0.6'              , 'autocorr'  , 'rho_spat' , 0.6  ; ...
    'AR(1) rho = 0.9'              , 'autocorr'  , 'rho_spat' , 0.9  ; ...
    'log-normal ROI noise'         , 'nonnormal' , ''         , NaN  ; ...
    'curvature c = 1'              , 'curvature' , 'curv'     , 1.0  ; ...
    'curvature c = 3'              , 'curvature' , 'curv'     , 3.0  ; ...
    'curvature c = 5'              , 'curvature' , 'curv'     , 5.0  };

nR = size(rows,1);
res = nan(nR,6);   % [power, ICC, R2m, rho, C, typeI]
B   = 1000;        % replicats bootstrap (identique a run_volet2_pipeline)

fprintf('=== Tests de robustesse, %d replicats, plan G=%d n=%d R=%d ===\n\n', ...
    NREP, P.G, P.n, P.R);

tt=tic;
for c=1:nR
    Pc = P;
    if ~isempty(rows{c,3}), Pc.(rows{c,3}) = rows{c,4}; end
    sc = rows{c,2};

    acc = nan(NREP,6);
    parfor mc=1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(sc, mc, c));   % bloc 1-19
        [X,Y] = gen_volet1(Pc, sc, rs);
        % bootstrap restreint a la paire de reference : la consistance de signe
        % n'est lue qu'au pic, comme dans run_volet2_pipeline, ce qui evite de
        % refaire B ajustements pour chacune des S*K paires.
        out   = paper_pipeline(X, Y, Pc, struct('bootstrap',true,'B',B, ...
                    'bootPairs',[s_star kRef]));

        g1 = (out.p_fdr < Pc.q) & (out.eta2 >= Pc.eta2_thr);

        acc(mc,:) = [ mean(g1(bandRef,kRef)) , ...   % puissance sur la bande
                      out.ICC(s_star,kRef)   , ...   % fiabilite au pic
                      out.R2m(s_star,kRef)   , ...   % NaN si verrou 2 non franchi
                      out.rho(s_star,kRef)   , ...   % Spearman au pic
                      out.C(s_star,kRef)     , ...   % NaN si verrou 3 non franchi
                      mean(g1(:,kNull))      ];      % type I sur le bruit ROI
    end
    res(c,:) = mean(acc,1,'omitnan');

    fprintf('%-24s Pow=%.3f  ICC=%.3f  R2=%.3f  rho=%.3f  C=%.3f  typeI=%.4f\n', ...
        rows{c,1}, res(c,1), res(c,2), res(c,3), res(c,4), res(c,5), res(c,6));
end
fprintf('\nTermine en %.1f min.\n', toc(tt)/60);

% --- lectures derivees, telles qu'annoncees dans le texte --------------------
iRef = 1;
fprintf('\n=== Lectures ===\n');
fprintf('Variation relative de R2 sous courbure (vs reference) :\n');
for c=6:8
    fprintf('  %-24s %+.1f %% (rho : %+.1f %%)\n', rows{c,1}, ...
        100*(res(c,3)-res(iRef,3))/res(iRef,3), ...
        100*(res(c,4)-res(iRef,4))/res(iRef,4));
end
fprintf('Erreur de type I sous autocorrelation : %.4f -> %.4f -> %.4f (ref %.4f)\n', ...
    res(2,6), res(3,6), res(4,6), res(iRef,6));
fprintf('Consistance bootstrap, minimum sur les huit lignes : %.3f\n', min(res(:,5)));

% --- sortie LaTeX, classe iopart --------------------------------------------
fid = fopen('stress_tests_table.tex','w');
fprintf(fid, '\\begin{tabular}{@{}lcccccc@{}}\n\\br\n');
fprintf(fid, ['Scenario & Power & Mean ICC & $R^2_{\\mathrm{marg}}$ & ' ...
              '$\\rho$ & $C$ & Type I error \\\\\n\\mr\n']);
for c=1:nR
    fprintf(fid, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.4f \\\\\n', ...
        rows{c,1}, res(c,1), res(c,2), res(c,3), res(c,4), res(c,5), res(c,6));
end
fprintf(fid, '\\br\n\\end{tabular}\n');
fclose(fid);
fprintf('\nTableau LaTeX ecrit dans stress_tests_table.tex\n');

save('ws_run_stress_tests.mat');
