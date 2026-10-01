% version v2, 10/09/2026 : le bruit porte sur le descripteur observe et non sur
%                          la reponse. La v1 ajoutait le bruit a Y, ce qui
%                          elargit le residu sans attenuer la pente : le
%                          descripteur restait le predicteur exact et le
%                          scenario ne testait pas ce qu'il annoncait. Les
%                          resultats produits par la v1 sont caducs.
% version v1, 10/09/2026 : z0 hors du cas exact
%RUN_LATENT  La correction de biais z0 discrimine-t-elle encore quand aucun
%   descripteur n'est la mesure exacte de ce qui produit la reponse.
%
%   LE PROBLEME. Dans le jeu de reference, la section 4.4 constate un contraste
%   net : z0 vaut 0.04 pour D1 a l'echelle generatrice et 1.04 pour D2 a la
%   meme echelle, alors que leurs R2 marginaux sont proches. La lecture
%   proposee est que z0 separe le descripteur qui porte la relation de celui
%   qui l'accompagne.
%
%   CE QUI EST SUSPECT. Le generateur construit Y directement sur
%   X(:,:,s_star,1). Il existe donc une paire pour laquelle la relation est
%   algebriquement exacte, sans terme d'erreur sur le predicteur. Le
%   reechantillonnage avec remise reduit la variation de X, ce qui attenue la
%   pente de toute paire ou X n'est qu'un proxy, mais laisse inchangee celle de
%   la paire exacte, X et Y se contractant ensemble. Le z0 nul de D1 pourrait
%   donc etre une propriete du generateur et non de l'indicateur.
%
%   Une surface reelle ne se comporte pas ainsi : aucun descripteur n'est la
%   cause du comportement fonctionnel, tous mesurent imparfaitement la
%   microgeometrie qui le produit. Si le contraste depend de l'exactitude, il
%   ne se retrouvera pas sur donnees mesurees et l'interpretation de la
%   section 4.4 doit etre restreinte.
%
%   PROTOCOLE. Le scenario 'latent' fixe la reponse sur le descripteur tel
%   qu'il est engendre, puis ajoute un bruit de mesure a tous les descripteurs
%   observes. Aucune paire n'est donc le predicteur exact, et toutes subissent
%   une attenuation au reechantillonnage. P.lat_noise regle l'ecart-type de ce
%   bruit, en fraction de l'ecart-type du descripteur porteur. A 0 le scenario
%   redevient identique a 'nominal' et sert de temoin.
%
%   Le bruit affecte aussi les deux premieres portes, ce qui est voulu : un
%   descripteur mesure avec erreur discrimine moins bien. Les colonnes de R2 et
%   de lambda du second tableau permettent de suivre cette degradation.
%
%   CE QUE LE BALAYAGE TRANCHE.
%     - Si z0 de la paire porteuse monte avec lat_noise et rejoint celui des
%       autres, le contraste est un artefact du cas exact.
%     - Si l'ecart entre porteuse et accompagnatrice se maintient, z0 ordonne
%       les descripteurs selon leur proximite a la relation, et la lecture de
%       la section 4.4 se generalise sous une forme attenuee.
%
%   Necessite : gen_volet1 v9, paper_pipeline, paper_params v9, paper_seed.
clear; clc;
pool = gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP  = 200;
LGRID = [0 0.25 0.50 0.75 1.00];   % bruit latent, en fraction de sd(X)
BOOT  = 1000;

P0 = paper_params();
ss = P0.s_star;                    % echelle generatrice
nL = numel(LGRID);

z0P = nan(nL,1); z0A = nan(nL,1);  % porteuse D1, accompagnatrice D2
r2P = z0P; r2A = z0P;
laP = z0P; laA = z0P;
ecart = z0P; fracSep = z0P;

fprintf('=== z0 hors du cas exact, %d replicats par point ===\n', NREP);
fprintf('Paire porteuse D1 a s = %d, accompagnatrice D2 a la meme echelle.\n\n', ss);
tt = tic;

for a = 1:nL
    P = P0; P.lat_noise = LGRID(a);
    zP = nan(NREP,1); zA = zP; rP = zP; rA = zP; lP = zP; lA = zP;
    parfor rep = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('nominal', rep, 270+a));
        [X,Y] = gen_volet1(P,'latent',rs);
        o = paper_pipeline(X,Y,P, struct('bootstrap',true,'B',BOOT, ...
                                         'bootPairs',[ss 1; ss 2]));
        zP(rep) = o.b1_z0(ss,1);  zA(rep) = o.b1_z0(ss,2);
        rP(rep) = o.R2m(ss,1);    rA(rep) = o.R2m(ss,2);
        lP(rep) = o.b1_ratio(ss,1); lA(rep) = o.b1_ratio(ss,2);
    end
    z0P(a)=median(zP,'omitnan'); z0A(a)=median(zA,'omitnan');
    r2P(a)=median(rP,'omitnan'); r2A(a)=median(rA,'omitnan');
    laP(a)=median(lP,'omitnan'); laA(a)=median(lA,'omitnan');
    ecart(a) = z0A(a) - z0P(a);
    % fraction de replicats ou l'ordre est respecte, porteuse sous
    % accompagnatrice : c'est la question pratique, pas l'ecart moyen
    v = ~isnan(zP) & ~isnan(zA);
    fracSep(a) = mean(zP(v) < zA(v));
    fprintf('lat_noise = %.2f termine (%.1f min)\n', LGRID(a), toc(tt)/60);
end

fprintf('\n--- Correction de biais z0 a l echelle generatrice ---\n');
fprintf('%10s %10s %10s %10s %14s\n', ...
    'lat_noise','z0 D1','z0 D2','ecart','ordre correct');
for a = 1:nL
    fprintf('%10.2f %10.3f %10.3f %10.3f %13.1f %%\n', ...
        LGRID(a), z0P(a), z0A(a), ecart(a), 100*fracSep(a));
end

fprintf('\n--- Contexte : R2 marginal et rapport lambda ---\n');
fprintf('%10s %10s %10s %10s %10s\n','lat_noise','R2 D1','R2 D2','lam D1','lam D2');
for a = 1:nL
    fprintf('%10.2f %10.3f %10.3f %10.2f %10.2f\n', ...
        LGRID(a), r2P(a), r2A(a), laP(a), laA(a));
end

fprintf(['\nLECTURE. La ligne lat_noise = 0 reproduit le cas exact du jeu de\n' ...
         'reference et sert de temoin. Si l ecart s effondre quand lat_noise\n' ...
         'augmente, le contraste de la section 4.4 tient a l exactitude du\n' ...
         'generateur. S il se maintient, meme reduit, z0 ordonne les paires\n' ...
         'selon leur proximite a la relation. La colonne ordre correct donne\n' ...
         'la part de replicats ou la porteuse a le z0 le plus bas : c est la\n' ...
         'question pratique, un ecart median favorable ne servant a rien s il\n' ...
         'ne tient pas replicat par replicat.\n']);

save('ws_run_latent.mat','LGRID','z0P','z0A','ecart','fracSep', ...
     'r2P','r2A','laP','laA','NREP','BOOT','P0');
fprintf('\nws_run_latent.mat ecrit. Total %.1f min.\n', toc(tt)/60);
