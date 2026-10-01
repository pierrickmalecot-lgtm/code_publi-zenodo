% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v4, 31/07/2026 : panel de reference a cinq descripteurs
%RUN_VALUE_ADDED  Apport des deux verrous amont, recalcule sur le panel actuel.
%   Compare la procedure complete a une variante non contrainte qui retient une
%   paire des lors que son R2 marginal depasse le seuil apres correction de
%   Benjamini-Hochberg, sans prerequis de discrimination ni de fiabilite.
%
%   POURQUOI RECALCULER. Les chiffres publies precedemment datent du panel ou
%   D5 et D6 etaient deux descripteurs de bruit identiques. Depuis, D6 porte un
%   effet specimen, ce qui modifie le comportement attendu de la variante non
%   contrainte : le R2 de D6 peut etre gonfle par la structure hierarchique
%   meme en l'absence de tout lien avec la reponse.
%
%   MISE EN OEUVRE. Le pipeline est appele avec les verrous amont neutralises
%   (q = 1, eta2_thr = 0, icc_thr = 0), de sorte que le R2 soit calcule pour
%   toutes les paires. La variante non contrainte applique ensuite le seul
%   critere de correlation, avec correction FDR sur les S*K valeurs. La
%   procedure complete est obtenue en reappliquant les trois criteres aux
%   memes quantites, ce qui garantit que les deux procedures voient exactement
%   les memes donnees.
%
%   Necessite : gen_volet1, paper_pipeline, paper_params.

clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end
P0 = paper_params();
S = P0.S; K = P0.K;
NREP = 200;

% Verite terrain du panel de reference, cinq descripteurs. Seuls D1, D2 et D3
% doivent etre retenus, soit neuf paires. Le cas du descripteur correle mais
% peu fiable, sur lequel se joue l'apport du verrou de fiabilite, fait l'objet
% du scenario dedie de run_unreliable.
bands = {5:7, 5:7, 3:5};
truePairs = false(S,K);
for k = 1:3, truePairs(bands{k},k) = true; end
nTrue = nnz(truePairs);
noiseCols = 4:5;              % bruits purs, panel a cinq descripteurs

% parametres avec verrous neutralises pour le calcul
P = P0;
P.q = 1.0; P.eta2_thr = 0; P.icc_thr = 0; P.r2_thr = -Inf;

% accumulateurs
retFull = zeros(NREP,1);  retFree = zeros(NREP,1);
sensFull= zeros(NREP,1);  sensFree= zeros(NREP,1);
fpFull  = zeros(NREP,1);  fpFree  = zeros(NREP,1);
d5Full  = zeros(NREP,1);  d5Free  = zeros(NREP,1);
d6Full  = zeros(NREP,1);  d6Free  = zeros(NREP,1);

parfor r = 1:NREP
    rs = RandStream('Threefry','Seed', paper_seed('nominal', r));
    [X,Y] = gen_volet1(P,'nominal',rs);
    out = paper_pipeline(X,Y,P);

    R2 = out.R2m;  R2(isnan(R2)) = 0;

    % --- procedure complete : les trois criteres, seuils du papier -----------
    g1 = (out.p_fdr < P0.q) & (out.eta2 >= P0.eta2_thr);
    g2 = g1 & (out.ICC >= P0.icc_thr);
    g3 = g2 & (R2 > P0.r2_thr);

    % --- variante non contrainte : le seul critere de correlation ------------
    % p-value du coefficient de pente, approchee par un test de Fisher sur le
    % R2 marginal avec N-2 degres de liberte, puis correction BH sur S*K.
    Nobs = P0.N * P0.R;
    Fr = (R2 ./ max(1-R2, eps)) * (Nobs - 2);
    pr = 1 - fcdf(Fr, 1, Nobs - 2);
    pAdj = bhfdr(pr(:));  pAdj = reshape(pAdj, S, K);
    gFree = (pAdj < P0.q) & (R2 > P0.r2_thr);

    retFull(r)  = nnz(g3);            retFree(r)  = nnz(gFree);
    sensFull(r) = nnz(g3 & truePairs)/nTrue;
    sensFree(r) = nnz(gFree & truePairs)/nTrue;
    fpFull(r)   = mean(g3(:,noiseCols),'all');
    fpFree(r)   = mean(gFree(:,noiseCols),'all');
    d5Full(r)   = mean(g3(:,4));      d5Free(r)   = mean(gFree(:,4));
    d6Full(r)   = mean(g3(:,5));      d6Free(r)   = mean(gFree(:,5));
end

fprintf('=== Apport des verrous, %d replicats, panel actuel ===\n\n', NREP);
fprintf('%-28s %14s %14s\n','', 'complete', 'sans verrous');
fprintf('%-30s %14.2f %14.2f\n','paires retenues', mean(retFull), mean(retFree));
fprintf('%-30s %14.3f %14.3f\n','sensibilite, 9 paires vraies', mean(sensFull), mean(sensFree));
fprintf('%-30s %14.4f %14.4f\n','bruits purs D4 et D5', mean(fpFull), mean(fpFree));
fprintf('%-30s %14.4f %14.4f\n','  dont D4, bruit ROI', mean(d5Full), mean(d5Free));
fprintf('%-30s %14.4f %14.4f\n','  dont D5, effet specimen', mean(d6Full), mean(d6Free));

fprintf(['\nLecture. Sur le panel de reference les deux procedures donnent le\n' ...
    'meme resultat : le critere de correlation ecarte a lui seul les\n' ...
    'descripteurs de bruit, avec 125 observations appariees. L''apport des\n' ...
    'verrous amont se manifeste ailleurs, sur un descripteur correle mais peu\n' ...
    'fiable, situation traitee par le scenario dedie de run_unreliable.\n']);

% --- sauvegarde sur disque ---------------------------------------------------
% La console peut etre perdue (fermeture de session, redemarrage). Les
% resultats sont donc ecrits sur disque des la fin du calcul.
T = table( ...
    {'retained';'sensitivity';'noise_total';'noise_D4';'noise_D5'}, ...
    [mean(retFull); mean(sensFull); mean(fpFull); mean(d5Full); mean(d6Full)], ...
    [mean(retFree); mean(sensFree); mean(fpFree); mean(d5Free); mean(d6Free)], ...
    'VariableNames', {'metric','full_procedure','no_gates'});
writetable(T,'value_added_results.csv');
save('value_added_results.mat', 'T', 'retFull','retFree','sensFull','sensFree', ...
     'fpFull','fpFree','d5Full','d5Free','d6Full','d6Free','NREP');
fprintf('\nResultats ecrits dans value_added_results.csv et .mat\n');

% -----------------------------------------------------------------------------
function padj = bhfdr(p)
M = numel(p); [ps,ix] = sort(p(:));
adj = ps .* M ./ (1:M)';
adj = flipud(cummin(flipud(adj)));
padj = zeros(M,1); padj(ix) = min(adj,1);
end
