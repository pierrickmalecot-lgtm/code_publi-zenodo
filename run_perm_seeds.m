% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/09/2026 : stabilite du niveau a n = 5, entre jeux de graines
%RUN_PERM_SEEDS  Le niveau des tests permutes a n = 5 est-il mesurable.
%
%   POURQUOI CE SCRIPT. Deux mesures independantes du taux d'erreur de type I
%   de la permutation studentisee sur D4, sous heteroscedasticite et a n = 5,
%   ont donne 0.0735 puis 0.0605. Chacune porte sur 2000 tests, donc un
%   ecart-type nominal de 0.0049 ; les deux valeurs different entre elles de
%   pres de deux ecarts-types. Un seul run ne permet donc pas de conclure, et
%   la premiere lecture, qui placait un seuil de validite vers n = 20, tirait
%   un seuil d'un point instable.
%
%   CE QUE MESURE CE SCRIPT. NRUN jeux de graines DISJOINTS, chacun de NREP
%   replicats. On regarde deux choses :
%     - la moyenne sur l'ensemble, qui vaut NRUN*NREP*S tests et resserre
%       l'incertitude d'un facteur racine de NRUN ;
%     - la DISPERSION entre runs, comparee a l'ecart-type binomial attendu. Si
%       la dispersion observee excede nettement l'attendu, c'est que la mesure
%       porte une variabilite supplementaire, et aucun run isole n'est
%       interpretable.
%
%   Seule la taille n = 5 est traitee : c'est la seule en cause, et c'est aussi
%   la seule accessible a un lot experimental courant en metrologie des
%   surfaces, ou G*n specimens au total est deja un plan lourd.
%
%   Necessite : gen_volet1, paper_pipeline v8, paper_params.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NRUN   = 4;        % jeux de graines disjoints
NREP   = 200;      % replicats par jeu
B_PERM = 10000;
SCEN   = 'heteroscedastic';
METH   = {'branching','perm','permwelch'};
MLBL   = {'logigramme papier','permutation (F)','permutation (Welch)'};
KNULL  = [4 5];

P = paper_params();   % n = 5, valeurs du papier
S = P.S; NM = numel(METH); NK = numel(KNULL);

alpha_run = zeros(NRUN, NM, NK);   % taux par run
fp_run    = zeros(NRUN, NM, NK);   % paires nulles au verrou 1, par run

fprintf('=== Stabilite du niveau a n = %d, %d jeux x %d replicats ===\n\n', ...
    P.n, NRUN, NREP);
tt = tic;

for ir = 1:NRUN
    aN = zeros(NREP,NM,NK); fN = aN;
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(SCEN, r, 120+ir));   % bloc 120-139
        [X,Y] = gen_volet1(P,SCEN,rs);
        aN_r = zeros(NM,NK); fN_r = aN_r;
        for im = 1:NM
            o = paper_pipeline(X,Y,P, struct('gate1Test',METH{im}, ...
                                             'B_perm',B_PERM));
            gate1 = (o.p_fdr < P.q) & (o.eta2 >= P.eta2_thr);
            for ik = 1:NK
                kk = KNULL(ik);
                aN_r(im,ik) = mean(o.p_raw(:,kk) < 0.05);
                fN_r(im,ik) = sum(gate1(:,kk));
            end
        end
        aN(r,:,:)=aN_r; fN(r,:,:)=fN_r;
    end
    alpha_run(ir,:,:) = squeeze(mean(aN,1));
    fp_run(ir,:,:)    = squeeze(mean(fN,1));
    fprintf('jeu %d termine, %.1f min ecoulees\n', ir, toc(tt)/60);
end
fprintf('\nTotal %.1f min.\n', toc(tt)/60);

% --- sorties ------------------------------------------------------------------
se_run  = sqrt(0.05*0.95/(S*NREP));            % attendu binomial, un run
se_pool = sqrt(0.05*0.95/(S*NREP*NRUN));       % attendu binomial, ensemble

for ik = 1:NK
    kk = KNULL(ik);
    fprintf('\n--- D%d, taux d erreur de type I par jeu de graines ---\n', kk);
    fprintf('Ecart-type binomial attendu : %.4f par jeu, %.4f sur l ensemble.\n', ...
        se_run, se_pool);
    fprintf('%-22s', 'methode');
    for ir = 1:NRUN, fprintf('   jeu %d', ir); end
    fprintf('    moyenne   ecart-type   observe/attendu   ecart au nominal\n');
    for im = 1:NM
        v = alpha_run(:,im,ik);
        fprintf('%-22s', MLBL{im});
        for ir = 1:NRUN, fprintf('  %6.4f', v(ir)); end
        z = (mean(v)-0.05)/se_pool;
        fprintf('   %8.4f   %10.4f   %15.2f   %+15.2f SE\n', ...
            mean(v), std(v), std(v)/se_run, z);
    end
end

fprintf('\n--- Paires nulles franchissant le verrou 1, moyenne sur %d echelles ---\n', S);
fprintf('%-22s   %12s   %12s\n','methode','D4','D5');
for im = 1:NM
    fprintf('%-22s   %12.3f   %12.3f\n', MLBL{im}, ...
        mean(fp_run(:,im,1)), mean(fp_run(:,im,2)));
end

fprintf('\nLECTURE. Une colonne observe/attendu proche de 1 signifie que la\n');
fprintf('dispersion entre jeux s explique par le seul alea binomial, et que la\n');
fprintf('moyenne sur l ensemble est interpretable. Nettement au-dessus de 1,\n');
fprintf('la mesure porte une variabilite supplementaire et aucun run isole ne\n');
fprintf('doit etre interprete.\n');

save('ws_run_perm_seeds.mat','alpha_run','fp_run','se_run','se_pool', ...
     'NRUN','NREP','B_PERM','SCEN','METH','MLBL','KNULL','P');
fprintf('\nws_run_perm_seeds.mat ecrit.\n');
