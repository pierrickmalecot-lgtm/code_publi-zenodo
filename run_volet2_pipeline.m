% version v10, 19/09/2026 : la coherence de signe C, saturee a 1 sur toute la
%                           grille, cede la place a la precision relative de la
%                           pente, qui varie avec le plan
% version v9, 04/09/2026 : ajout de rho_nonmono, absent car Pc est construit a la main
% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v3, 01/08/2026 : dimensions lues dans paper_params, taux d erreur mesure sur le bruit ROI
%RUN_VOLET2_PIPELINE  Etude de dimensionnement (5.9) avec le PIPELINE REEL.
%   Remplace run_volet2.m, qui reimplementait une version allegee du verrou 1
%   (F poole sans FDR, sans branches robustes, sans chainage des verrous) et
%   son propre generateur. Ici, une seule implementation fait foi :
%   gen_volet1 pour la generation, paper_pipeline pour l'analyse.
%
%   Grille : G x n x R x niveau de bruit (81 configurations).
%   Metriques, toutes mesurees SUR LA CHAINE COMPLETE :
%     power  : frequence de passage du verrou 1 (p_FDR < q ET eta2 >= seuil)
%              pour D1 sur sa bande fonctionnelle 5-7
%     ICC    : ICC moyen de D1 a l'echelle de reference s*
%     R2     : R2 marginal de D1 a s* (NaN si la paire n'atteint pas l'etage)
%     rp     : precision relative de la pente a la paire de reference, soit
%              l'ecart-type bootstrap de beta1 rapporte a sa valeur absolue
%              (NaN si la paire n'atteint pas l'etage de correlation). Elle
%              remplace la coherence de signe C, qui valait 1 sur les 81
%              configurations et n'apprenait donc rien ici ; elle est en outre
%              la quantite que la section sur l'ecriture d'une specification
%              demande de reporter.
%     FP     : frequence de passage du verrou 1 pour D5 (bruit ROI, null propre)
%
%   COUT : le pipeline traite les S*K = 60 paires a chaque replicat, contre une
%   seule echelle dans l'ancienne version. Le bootstrap est donc restreint a la
%   paire de reference (opts.bootPairs), ce qui ramene le surcout a un facteur
%   1,5 a 3 environ. Reduire nMC ou B si necessaire.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end

nMC = 200;          % reduit par rapport a 500 : voir note de cout
B   = 200;          % replicats bootstrap, paire de reference seulement

G_levels=[3 5 7]; n_levels=[3 5 10]; R_levels=[3 5 10];
noise_levels={'low','medium','high'};
% sb : ecart-type specimen | sw : ROI | sf : bruit fonctionnel
% amp : facteur d'amplitude du signal condition | b1 : pente vraie
np.low   =struct('sb',0.10,'sw',0.08,'sf',0.25,'amp',0.70,'b1',2.2);
np.medium=struct('sb',0.18,'sw',0.15,'sf',0.40,'amp',0.55,'b1',1.8);
np.high  =struct('sb',0.30,'sw',0.25,'sf',0.65,'amp',0.45,'b1',1.4);

% Dimensions lues dans paper_params, de sorte qu'un changement de panel n'ait
% pas a etre reporte ici. Le null de reference est D4, bruit ROI seul : D5
% porte un effet specimen et franchit le premier verrou par pseudo-replication,
% ce n'est donc pas un null propre.
Pref = paper_params();
S = Pref.S; K = Pref.K; s_star = Pref.s_star;
kRef = 1; bandRef = 5:7; kNull = 4;

cfgs={}; ci=0;
for iG=1:3, for iN=1:3, for iR=1:3, for iZ=1:3
    ci=ci+1; cfgs{ci}=struct('G',G_levels(iG),'n',n_levels(iN), ...
        'R',R_levels(iR),'nz',noise_levels{iZ},'p',np.(noise_levels{iZ}));
end,end,end,end
nC=numel(cfgs);
res=zeros(nC,5);   % pow, icc, r2, rp, fp
tt=tic;

for c=1:nC
    cfg=cfgs{c}; q=cfg.p;
    % --- parametres du modele pour cette configuration -----------------------
    Pc = struct();
    Pc.G=cfg.G; Pc.n=cfg.n; Pc.R=cfg.R; Pc.S=S; Pc.K=K; Pc.N=cfg.G*cfg.n;
    Pc.condLevels   = q.amp * linspace(0.2,1.0,cfg.G);   % amplitude du signal
    Pc.sigma_between= q.sb;
    Pc.sigma_within = q.sw;
    Pc.sigma_func   = q.sf;
    Pc.sigma_u      = 0.20;
    Pc.beta0=2.0; Pc.beta1=q.b1; Pc.s_star=s_star;
    % Reglages du terme local. Ils doivent figurer ici car cette etude
    % construit ses parametres condition par condition au lieu d'appeler
    % paper_params, et le generateur les reclame. Ils restent sans effet, le
    % terme local n'existant que sous le scenario 'unreliable'.
    Pc.sigma_z = Pref.sigma_z;
    Pc.gamma_z = Pref.gamma_z;
    % Meme raison pour le seuil du drapeau de non-monotonie, que le pipeline
    % lit systematiquement. Il ne change aucune decision : le drapeau ne
    % s'applique qu'aux paires deja rejetees par la porte 3.
    Pc.rho_nonmono = Pref.rho_nonmono;
    Pc.q=0.05; Pc.eta2_thr=0.06; Pc.icc_thr=0.60; Pc.r2_thr=0.30; Pc.C_thr=0.90;

    acc=nan(nMC,5);
    parfor mc=1:nMC
        rs=RandStream('Threefry','Seed', paper_seed('nominal', mc, 300+c));   % bloc 300-399
        [X,Y]=gen_volet1(Pc,'nominal',rs);
        out=paper_pipeline(X,Y,Pc, struct('bootstrap',true,'B',B, ...
            'bootPairs',[s_star kRef]));

        g1 = (out.p_fdr < Pc.q) & (out.eta2 >= Pc.eta2_thr);

        pow  = mean(g1(bandRef,kRef));      % puissance sur la bande vraie
        icc  = out.ICC(s_star,kRef);
        r2   = out.R2m(s_star,kRef);        % NaN si verrou 2 non franchi
        rp   = out.b1_rp(s_star,kRef);      % NaN si verrou 3 non franchi
        fp   = mean(g1(:,kNull));           % faux positifs sur le bruit ROI

        acc(mc,:)=[pow, icc, r2, rp, fp];
    end
    res(c,:)=mean(acc,1,'omitnan');
    fprintf(['%3d/%d [G=%d n=%2d R=%2d %6s] Pow=%.2f ICC=%.2f R2=%.2f ' ...
        'rp=%.3f FP=%.3f (ETA %.0f s)\n'], ...
        c,nC,cfg.G,cfg.n,cfg.R,cfg.nz,res(c,:),toc(tt)/c*(nC-c));
end

% La table est construite en une seule fois. L'assemblage precedent, colonne par
% colonne et ligne par ligne, faisait grandir la table a chaque nouvelle ligne
% avant que toutes ses colonnes soient renseignees, ce qui declenchait un
% avertissement par configuration, quatre-vingts au total dans le journal. Les
% valeurs produites etaient correctes, toutes les colonnes etant assignees dans
% la meme iteration.
G_ = zeros(nC,1); n_ = G_; R_ = G_; nz_ = strings(nC,1);
for c = 1:nC
    G_(c) = cfgs{c}.G; n_(c) = cfgs{c}.n; R_(c) = cfgs{c}.R;
    nz_(c) = string(cfgs{c}.nz);
end
T = table(G_, n_, R_, nz_, res(:,1), res(:,2), res(:,3), res(:,4), res(:,5), ...
    'VariableNames', {'G','n','R','noise','power_F','mean_ICC', ...
                      'R2_recovery','rel_precision','false_positive_rate'});
save('volet2_pipeline_results.mat','T');
writetable(T,'volet2_results_paper.csv');   % meme nom : figures inchangees
fprintf(['Termine en %.1f min. CSV ecrit sous volet2_results_paper.csv,\n' ...
    'compatible avec make_sensitivity_figures.\n'], toc(tt)/60);
