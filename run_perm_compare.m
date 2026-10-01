% version v9, 27/09/2026 : controle ROI en Welch (permobswelch) ; variantes F
%                         retirees, le papier ne presente que Welch
% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v3, 02/09/2026 : ajout de la variante permwelch
% version v2, 02/09/2026 : deux nulls (D4, D5), verrou 1 ventile par descripteur
% version v1, 02/09/2026 : logigramme conditionnel contre permutation
%RUN_PERM_COMPARE  L'ANOVA par permutation inconditionnelle peut-elle remplacer
%   le branchement Shapiro-Wilk / Brown-Forsythe / Welch / Kruskal-Wallis.
%
%   L'ALTERNATIVE. Le verrou 1 du papier route la p-value selon deux tests
%   preliminaires : normalite par Shapiro-Wilk, homogeneite par Brown-Forsythe.
%   L'objection faite a ce schema est double. D'une part, a n petit ces deux
%   tests ont une puissance faible, donc le routage est lui-meme bruite et le
%   resultat final depend de deux decisions peu fiables. D'autre part, un test
%   conditionnel a un test preliminaire n'a pas le niveau nominal de l'un ni de
%   l'autre. L'ANOVA par permutation ne demande ni normalite ni homogeneite et
%   fournit une p-value unique, sans branchement.
%
%   UNITE DE PERMUTATION. Point decisif. Sous H0, ce sont les SPECIMENS qui
%   sont echangeables entre conditions : les R regions mesurees sur un meme
%   specimen partagent son effet aleatoire, donc ne sont pas echangeables une a
%   une. Ce script compare les deux mises en oeuvre, 'perm' au niveau specimen
%   et 'permobswelch' au niveau observation, pour mesurer ce que coute l'erreur.
%
%   TROIS METHODES COMPAREES (les variantes en F de la v8 sont retirees,
%   le papier ne presente que Welch).
%     branching     logigramme du papier. Le resultat de Brown-Forsythe n'est
%                   consulte que si la normalite est acceptee, donc des donnees
%                   a la fois non normales et heteroscedastiques partent en
%                   Kruskal-Wallis, qui ne corrige pas les variances inegales.
%     permwelch     permutation des specimens, statistique de Welch.
%                   Studentiser la statistique rend le test valide sous
%                   variances inegales : c'est la seule variante qui remplace
%                   les DEUX tests preliminaires.
%     permobswelch  permutation des observations, statistique de Welch. Unite
%                   de permutation incorrecte, fournie comme controle : a
%                   statistique egale, seule l'unite change.
%
%   DEUX NULLS, DEUX QUESTIONS DIFFERENTES.
%     ATTENTION. Le scenario 'null' de gen_volet1 supprime le lien X vers Y :
%     il est le null du verrou 3, PAS celui du verrou 1. L'effet de condition
%     sur X y subsiste, donc D1, D2 et D3 y discriminent legitimement les
%     conditions. Les nulls du verrou 1 sont des COLONNES du panel, presentes
%     dans tous les scenarios.
%
%     D4, bruit de region pur, sans effet condition ni effet specimen.
%         Null de reference classique : mesure le niveau nominal du test.
%         Aucune structure de grappe n'y existe, donc ce descripteur ne peut
%         PAS reveler un probleme d'echangeabilite.
%
%     D5, bruit specimen plus bruit region, sans effet condition.
%         C'est le null qui teste l'unite de permutation. L'effet specimen
%         gonfle MSbetween dans le pooling n*R, donc le F standard y depasse 1
%         sous H0 : c'est la pseudo-replication au sens de Hurlbert. Une
%         permutation au niveau SPECIMEN deplace les specimens en bloc, sa
%         distribution de reference contient donc cette variance inter-specimens
%         et le test en tient compte. Une permutation au niveau OBSERVATION
%         casse les grappes et ne le fait pas. L'ecart entre les deux sur D5
%         mesure directement ce que coute l'erreur d'unite.
%
%   VENTILATION. Les compteurs du verrou 1 sont donnes par descripteur, pour
%   distinguer les paires legitimement retenues (D1 a D3 dans leur bande) de
%   celles qui ne le sont que par pseudo-replication (D5, et hors bande).
%
%     null             absence de lien X vers Y : verifie qu'aucune paire ne
%                      ressort en sortie de pipeline.
%     nominal          puissance et accord des decisions avec le branchement.
%     nonnormal        cas ou le branchement route vers Kruskal-Wallis.
%     heteroscedastic  cas ou le branchement route vers Welch.
%
%   Necessite : gen_volet1, paper_pipeline v11 (variante permobswelch),
%   paper_params, paper_seed.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP   = 200;
B_PERM = 10000;   % resolution 1e-4, sous le plus petit seuil BH q/M = 1e-3
KNULL  = [4 5];   % colonnes nulles : bruit ROI pur, puis bruit specimen + ROI
SCEN   = {'null','nominal','nonnormal','heteroscedastic'};
METH   = {'branching','permwelch','permobswelch'};
MLBL   = {'logigramme papier','permutation specimen (Welch)', ...
          'permutation observation (Welch)'};

P = paper_params();
S = P.S; K = P.K;
NS = numel(SCEN); NM = numel(METH); NK = numel(KNULL);

alpha_null = zeros(NS, NM, NK);   % p_raw < 0.05 sur chaque colonne nulle
fp_null    = zeros(NS, NM, NK);   % paires nulles franchissant le verrou 1
any_null   = zeros(NS, NM, NK);   % replicats avec au moins une
g1_desc    = zeros(NS, NM, K);    % verrou 1 ventile par descripteur
n_gate1    = zeros(NS, NM);       % verrou 1, toutes colonnes
n_gate3    = zeros(NS, NM);       % paires retenues en sortie
agree      = zeros(NS, NM);       % accord verrou 1 avec branching

fprintf('=== Branchement contre permutation, %d replicats, B_perm = %d ===\n\n', ...
    NREP, B_PERM);
tt = tic;

for is = 1:NS
    scen = SCEN{is};
    aN = zeros(NREP,NM,NK); fN = aN; yN = aN;
    gD = zeros(NREP,NM,K);
    g1 = zeros(NREP,NM); ag = g1; g3 = g1;

    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(scen, r));
        [X,Y] = gen_volet1(P,scen,rs);
        aN_r = zeros(NM,NK); fN_r = aN_r; yN_r = aN_r;
        gD_r = zeros(NM,K);
        g1_r = zeros(1,NM); ag_r = g1_r; g3_r = g1_r;
        ref = [];
        for im = 1:NM
            o = paper_pipeline(X,Y,P, struct('gate1Test',METH{im}, ...
                                             'B_perm',B_PERM));
            gate1 = (o.p_fdr < P.q) & (o.eta2 >= P.eta2_thr);
            if im == 1, ref = gate1; end
            for ik = 1:NK
                kk = KNULL(ik);
                aN_r(im,ik) = mean(o.p_raw(:,kk) < 0.05);
                fN_r(im,ik) = sum(gate1(:,kk));
                yN_r(im,ik) = double(any(gate1(:,kk)));
            end
            gD_r(im,:) = sum(gate1,1);
            g1_r(im)   = sum(gate1(:));
            ag_r(im)   = mean(gate1(:) == ref(:));
            g3_r(im)   = sum(o.decision(:) == 3);
        end
        aN(r,:,:)=aN_r; fN(r,:,:)=fN_r; yN(r,:,:)=yN_r; gD(r,:,:)=gD_r;
        g1(r,:)=g1_r; ag(r,:)=ag_r; g3(r,:)=g3_r;
    end

    alpha_null(is,:,:) = squeeze(mean(aN,1));
    fp_null(is,:,:)    = squeeze(mean(fN,1));
    any_null(is,:,:)   = squeeze(mean(yN,1));
    g1_desc(is,:,:)    = squeeze(mean(gD,1));
    n_gate1(is,:)      = mean(g1,1);
    n_gate3(is,:)      = mean(g3,1);
    agree(is,:)        = mean(ag,1);
    fprintf('%-16s termine (%.1f min ecoulees)\n', scen, toc(tt)/60);
end
fprintf('\nTotal %.1f min.\n', toc(tt)/60);

% --- sorties ------------------------------------------------------------------
for ik = 1:NK
    kk = KNULL(ik);
    if kk == 4
        titre = 'D4, bruit de region pur : niveau nominal du test';
    else
        titre = 'D5, bruit specimen + region : test de l unite de permutation';
    end
    fprintf('\n--- %s ---\n', titre);
    fprintf('Niveau nominal : 0.05 pour p_raw brut, %.2f pour le verrou 1 apres FDR.\n', P.q);
    fprintf('%-16s %-26s  p_raw<0.05   paires verrou 1   replicats avec au moins une\n', ...
        'scenario','methode');
    for is = 1:NS
        for im = 1:NM
            fprintf('%-16s %-26s  %9.4f   %15.3f   %26.3f\n', SCEN{is}, MLBL{im}, ...
                alpha_null(is,im,ik), fp_null(is,im,ik), any_null(is,im,ik));
        end
    end
end

fprintf('\n--- Verrou 1 ventile par descripteur (paires sur %d echelles) ---\n', S);
for is = 1:NS
    fprintf('\n%s\n', SCEN{is});
    fprintf('%-26s', 'methode');
    for k = 1:K, fprintf('     D%d', k); end
    fprintf('   total  retenues  accord\n');
    for im = 1:NM
        fprintf('%-26s', MLBL{im});
        for k = 1:K, fprintf('  %5.2f', g1_desc(is,im,k)); end
        fprintf('  %6.2f  %8.2f  %6.3f\n', n_gate1(is,im), n_gate3(is,im), agree(is,im));
    end
end

% table longue : une ligne par (scenario, methode)
sc = repelem(SCEN(:),NM,1); me = repmat(MLBL(:),NS,1);
T = table(sc, me, reshape(alpha_null(:,:,1)',[],1), reshape(fp_null(:,:,1)',[],1), ...
    reshape(alpha_null(:,:,2)',[],1), reshape(fp_null(:,:,2)',[],1), ...
    reshape(n_gate1',[],1), reshape(n_gate3',[],1), reshape(agree',[],1), ...
    'VariableNames',{'scenario','method','alpha_D4','fp_D4','alpha_D5','fp_D5', ...
                     'n_gate1','n_gate3','agreement'});
for k = 1:K
    T.(sprintf('gate1_D%d',k)) = reshape(g1_desc(:,:,k)',[],1);
end
writetable(T,'perm_compare_results.csv');
save('ws_run_perm_compare.mat','alpha_null','fp_null','any_null','g1_desc', ...
     'n_gate1','n_gate3','agree','SCEN','METH','MLBL','KNULL','NREP','B_PERM','P','T');
fprintf('\nperm_compare_results.csv et ws_run_perm_compare.mat ecrits.\n');
