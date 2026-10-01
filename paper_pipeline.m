% version v11, 27/09/2026 : variante permobswelch (observations, Welch), controle
%                          de l'unite de permutation a statistique egale
% version v10, 27/09/2026 : correction n/(n-1) de la variance bootstrap
%                          (opts.nn1, active par defaut) ; intervalle et
%                          ecart-type sans correction conserves (b1_*_raw)
% version v9, 04/09/2026 : drapeau de non-monotonie et reajustement quadratique
% version v8, 02/09/2026 : variante permwelch, statistique de Welch permutee
% version v7, 02/09/2026 : correction de l ordre de mise a plat dans la variante
%                          de permutation par observation (perm_pvalue)
% version v6, 02/09/2026 : option gate1Test, ANOVA par permutation
% version v5, 02/09/2026 : intervalle BCa en remplacement du percentile,
%                          ecart-type parametrique de la pente et rapport
%                          bootstrap / parametrique
% version v4, 02/09/2026 : boot_C renvoie le vecteur des pentes, denominateur
%                          corrige des ajustements non aboutis, standardisation
%                          bootstrap fixee sur l'echantillon complet
% version v3, 31/07/2026 : options diagMode et bootPairs, ventilation des branches
function out = paper_pipeline(X, Y, P, opts)
%PAPER_PIPELINE  Pipeline STRICTEMENT conforme au papier (sections 3.6-3.9).
%   - Verrous 1-2 : POOLING n*R par condition (X_{g,j}, j=1..m=n*R)
%     * MS/F : eqs msbetween/mswithin/F ; eta2/omega2 : eq effectsizes
%     * Brown-Forsythe -> Welch ; (Shapiro-Wilk si swtest dispo)
%     * BH-FDR sur S*K (eq 25) ; gate1 : p_fdr<q & eta2>=thr
%     * ICC = (F-1)/(F+m-1), m=n*R (eq icc_f) ; gate2 : ICC>=thr
%   - Verrou 3 : LME Y~X+(1|specimen), R2 MARGINAL (eq 18) ; gate3 : R2>thr
%   - Spearman (eq spearman) sur les N_obs observations appariees
%   - opts.bootstrap : C_{s,k}, bootstrap specimen STRATIFIE par condition
%     Sorties du bootstrap : C (coherence de signe), Cboot (vecteur des pentes),
%     b1_lo / b1_hi (intervalle BCa), b1_pct_lo / b1_pct_hi (percentile brut,
%     conserve pour comparaison), b1_z0 et b1_accel (parametres BCa),
%     b1_se, b1_rp (precision relative), b1_ratio (ecart-type bootstrap rapporte
%     a l'ecart-type parametrique), Cinfo (comptes d'ajustements).
%   - out.b1_se_par : ecart-type de la pente rendu par fitlme, calcule a
%     echantillon de specimens fixe. Disponible des le verrou 3, sans bootstrap.
%   - opts.nn1 (defaut true) : la dispersion des pentes bootstrap est
%     multipliee par sqrt(n/(n-1)) autour de l'estimation ponctuelle avant le
%     calcul de C, de l'intervalle BCa, de b1_se, b1_rp et b1_ratio. Les champs
%     b1_lo_raw, b1_hi_raw, b1_pct_lo_raw, b1_pct_hi_raw et b1_se_raw gardent
%     les valeurs sans correction, pour comparaison.
%   - opts.loso : drapeau fragile (3.9)
if nargin<4, opts=struct; end
def = {'bootstrap',false;'B',1000;'loso',false;'diagMode','residual';'bootPairs',[];'bootAlpha',0.05; ...
       'gate1Test','branching';'B_perm',10000;'quadCheck',false;'nn1',true};
for d=1:size(def,1), if ~isfield(opts,def{d,1}), opts.(def{d,1})=def{d,2}; end, end
[N,R,S,K] = size(X); G=P.G; n=P.n; m=n*R;
sampleCondition = repelem(1:G, n);

F=nan(S,K); p_raw=F; eta2=F; omega2=F; ICC=F; R2m=F; rho=F; b1hat=F; Ftest=F;
b1se_par=nan(S,K);
branch=zeros(S,K);   % 0=F standard, 1=Welch, 2=Kruskal-Wallis, 3=permutation
                     % specimen (F), 4=permutation observation (F),
                     % 5=permutation specimen (Welch)
                     % 6=permutation observation (Welch)

% Permutations tirees une seule fois et partagees par toutes les paires (s,k).
% Partager le meme jeu preserve la structure de dependance entre les S*K tests,
% ce qui est coherent avec la correction FDR appliquee ensuite.
PI = [];
if ~strcmpi(opts.gate1Test,'branching')
    U = G*n;                                  % unites = specimens
    if any(strcmpi(opts.gate1Test,{'permobs','permobswelch'})), U = G*m; end
    [~,PI] = sort(rand(opts.B_perm, U), 2);
end

for s=1:S, for k=1:K
    % pooling X(g,j)
    Xp = zeros(G,m);
    for g=1:G
        Xp(g,:) = reshape(X(sampleCondition==g,:,s,k), 1, m);
    end
    [F(s,k), p_raw(s,k), eta2(s,k), omega2(s,k), Ftest(s,k), branch(s,k)] = ...
        anova_gate_paper(Xp, n, R, opts.diagMode, opts.gate1Test, PI);
    % ICC depuis le F STANDARD uniquement (la formule (F-1)/(F+m-1) n'est pas
    % valide pour le F de Welch ; la branche robuste ne sert qu'a la p-value)
    ICC(s,k) = max((F(s,k)-1)/(F(s,k)+m-1), 0);
    % Spearman sur observations appariees (mid-ranks)
    rho(s,k) = corr(reshape(X(:,:,s,k),[],1), Y(:), 'Type','Spearman');
end, end

p_fdr = reshape(bh_fdr(p_raw(:)), S, K);
gate1 = (p_fdr < P.q) & (eta2 >= P.eta2_thr);
gate2 = gate1 & (ICC >= P.icc_thr);

for s=1:S, for k=1:K
    if ~gate2(s,k), continue; end
    [R2m(s,k), b1hat(s,k), b1se_par(s,k)] = lme_r2_marginal(X(:,:,s,k), Y);
end, end
gate3 = gate2 & (R2m > P.r2_thr);

% --- drapeau de non-monotonie (inspection, pas une porte) ---------------------
% Une paire qui a franchi les portes 1 et 2 discrimine les conditions et est
% mesuree de facon reproductible : il y a bien une structure. Si son R2
% marginal est neanmoins sous le seuil ET que son Spearman est faible en
% valeur absolue, ce n'est pas l'absence de relation qui est etablie mais
% l'inadequation de la forme ajustee. Une relation a optimum interieur, en U ou
% en cloche, annule simultanement la pente lineaire et la correlation de rang.
% Sans ce signalement la paire est rejetee a la porte 3 sans que rien
% n'indique pourquoi, et l'extension quadratique de la section 3.9.6 n'est
% jamais declenchee.
%
% Ce drapeau ne modifie aucune decision. Il ne s'applique qu'aux paires deja
% rejetees par la porte 3, il ne peut donc ni en retenir une de plus ni en
% exclure une.
out_nonmono = gate2 & ~gate3 & (abs(rho) < P.rho_nonmono);

out = struct('F',F,'Ftest',Ftest,'p_raw',p_raw,'p_fdr',p_fdr,'eta2',eta2,...
             'omega2',omega2,'ICC',ICC,'R2m',R2m,'rho',rho,'beta1',b1hat);
out.b1_se_par = b1se_par;   % ecart-type parametrique, echantillon fixe
out.nonmono   = out_nonmono; % relation probablement non monotone (inspection)

% --- reajustement quadratique sur les paires signalees ------------------------
% opts.quadCheck = false par defaut : aucun resultat existant n'est modifie.
% Actif, le bloc ajuste Y ~ X + X^2 + (1|specimen) sur les seules paires
% portant le drapeau, et rapporte le R2 marginal du modele etendu ainsi que
% l'ecart d'AIC par rapport au modele lineaire. Un dAIC positif indique que le
% terme quadratique merite les deux degres de liberte qu'il coute, selon la
% regle de selection de la section 3.9.6. Le cout est d'un ajustement par
% paire signalee, aucune si le drapeau ne se leve pas.
if opts.quadCheck
    out.R2m_quad = nan(S,K);
    out.dAIC     = nan(S,K);
    for s=1:S, for k=1:K
        if out_nonmono(s,k)
            [out.R2m_quad(s,k), out.dAIC(s,k)] = lme_quad(X(:,:,s,k), Y);
        end
    end, end
end
out.branch = branch;   % 0=F standard, 1=Welch, 2=Kruskal-Wallis
% out.F = F standard (MSb/MSw), utilise pour ICC/effect sizes et les figures ;
% out.Ftest = statistique de la branche effectivement testee (F, Welch, KW).
out.decision = double(gate1)+double(gate2)+double(gate3);

if opts.bootstrap
    out.C     = nan(S,K);   % coherence de signe
    out.b1_lo = nan(S,K);   % borne basse BCa
    out.b1_hi = nan(S,K);   % borne haute BCa
    out.b1_pct_lo = nan(S,K);  % percentile brut, conserve pour comparaison
    out.b1_pct_hi = nan(S,K);
    out.b1_z0    = nan(S,K);   % correction de biais
    out.b1_accel = nan(S,K);   % acceleration
    out.b1_se = nan(S,K);   % ecart-type bootstrap de la pente
    out.b1_rp = nan(S,K);   % precision relative : b1_se / |b1hat|
    out.b1_ratio = nan(S,K);% b1_se / b1_se_par
    out.Cboot = cell(S,K);  % vecteur des B pentes, pour recalcul sans relance
    out.Cinfo = cell(S,K);  % comptes d'ajustements aboutis / non aboutis
    out.b1_method = cell(S,K);  % 'BCa' ou 'percentile' selon le domaine de validite
    out.b1_lo_raw = nan(S,K);   % memes grandeurs sans la correction n/(n-1)
    out.b1_hi_raw = nan(S,K);
    out.b1_pct_lo_raw = nan(S,K);
    out.b1_pct_hi_raw = nan(S,K);
    out.b1_se_raw = nan(S,K);
    out.b1_nn1 = opts.nn1;      % correction appliquee ou non
    % opts.bootPairs = [] -> toutes les paires retenues (defaut, conforme au
    % papier). Sinon matrice [p x 2] de couples (s,k) : restreint le bootstrap
    % a ces paires, ce qui divise le cout dans les etudes Monte Carlo ou seule
    % la paire de reference est suivie.
    if isempty(opts.bootPairs)
        sel = gate3;
    else
        sel = false(S,K);
        for p = 1:size(opts.bootPairs,1)
            sel(opts.bootPairs(p,1), opts.bootPairs(p,2)) = true;
        end
        sel = sel & gate3;
    end
    for s=1:S, for k=1:K
        if sel(s,k)
            [out.C(s,k), bb, out.Cinfo{s,k}] = ...
                boot_C(X(:,:,s,k), Y, b1hat(s,k), P, opts.B);
            out.Cboot{s,k} = bb;        % pentes brutes, sans correction
            v = bb(~isnan(bb));
            if numel(v) >= 2
                % VARIANCE n/(n-1). Le tirage avec remise de n specimens par
                % condition sous-estime la dispersion de la pente (Field et
                % Welsh 2007). Chaque pente est ecartee de l'estimation
                % ponctuelle du facteur sqrt(n/(n-1)) : la variance est
                % multipliee par n/(n-1) et la position des pentes par rapport
                % a l'estimation, donc z0, est inchangee. Le facteur tend vers
                % 1 quand n grandit (1.25 a n = 5, 1.11 a n = 10).
                if opts.nn1
                    fk = sqrt(n/(n-1));
                else
                    fk = 1;
                end
                vc = b1hat(s,k) + fk*(v - b1hat(s,k));
                ci = bca_ci(vc, b1hat(s,k), X(:,:,s,k), Y, opts.bootAlpha);
                % C recalcule sur le nuage corrige : l'elargissement peut faire
                % passer une pente proche de zero de l'autre cote.
                out.C(s,k)         = sum(sign(vc) == sign(b1hat(s,k))) / numel(vc);
                out.b1_lo(s,k)     = ci.lo;
                out.b1_hi(s,k)     = ci.hi;
                out.b1_pct_lo(s,k) = ci.pct_lo;
                out.b1_pct_hi(s,k) = ci.pct_hi;
                out.b1_z0(s,k)     = ci.z0;
                out.b1_accel(s,k)  = ci.accel;
                out.b1_se(s,k)     = std(vc);
                out.b1_rp(s,k)     = std(vc)/abs(b1hat(s,k));
                out.b1_ratio(s,k)  = std(vc)/b1se_par(s,k);
                out.b1_method{s,k} = ci.method;
                % valeurs sans correction, pour comparaison
                if fk == 1
                    ci0 = ci;
                else
                    ci0 = bca_ci(v, b1hat(s,k), X(:,:,s,k), Y, opts.bootAlpha);
                end
                out.b1_lo_raw(s,k)     = ci0.lo;
                out.b1_hi_raw(s,k)     = ci0.hi;
                out.b1_pct_lo_raw(s,k) = ci0.pct_lo;
                out.b1_pct_hi_raw(s,k) = ci0.pct_hi;
                out.b1_se_raw(s,k)     = std(v);
            end
        end
    end, end
end
if opts.loso, out.fragile = loso_flags(X, Y, out, P); end
end

% -----------------------------------------------------------------------------
function [Fstd,p,eta2,omega2,Fv,branch] = anova_gate_paper(Xp, n, R, diagMode, gate1Test, PI)
% diagMode : 'residual' (defaut) teste les hypotheses sur les residus intra-
%            specimen w = X - moyenne_specimen (coherent avec la structure
%            d'erreur du modele mixte, verrou 3) ;
%            'pooled' teste sur X poole brut (variante conservatrice : l'effet
%            specimen a_i est lu comme une violation -> route plus souvent en
%            robuste). Ne change QUE le routage F/Welch/KW, pas les statistiques
%            de discrimination (F, eta2, ICC) ni les verrous aval.
%
% gate1Test : 'branching' (defaut) applique le logigramme du papier, soit
%             Shapiro-Wilk puis Brown-Forsythe, puis F standard, Welch ou
%             Kruskal-Wallis selon le resultat. Le resultat de Brown-Forsythe
%             n'est consulte que si la normalite est acceptee : des donnees a
%             la fois non normales et heteroscedastiques partent en
%             Kruskal-Wallis, qui ne corrige pas l'inegalite des variances.
%
%             Les trois valeurs suivantes remplacent tout le branchement par un
%             test unique, applique sans routage conditionnel. Elles repondent
%             a l'objection selon laquelle un test choisi en fonction de deux
%             tests preliminaires peu puissants n'a le niveau nominal ni de
%             l'un ni de l'autre.
%
%             'perm'      permutation des specimens, statistique F ordinaire.
%             'permobs'   permutation des observations, statistique F. Variante
%                         fournie pour comparaison seulement : les R regions
%                         d'un meme specimen partagent son effet aleatoire et ne
%                         sont pas echangeables separement.
%             'permwelch' permutation des specimens, statistique de Welch.
%             'permobswelch' permutation des observations, statistique de
%                         Welch. Controle de permwelch a statistique egale :
%                         seule l'unite de permutation change.
%
%   POURQUOI LA STATISTIQUE DE WELCH. Un test par permutation teste
%   l'echangeabilite, pas l'egalite des moyennes. Sous variances inegales, les
%   observations n'ont pas la meme loi d'une condition a l'autre : elles ne sont
%   pas echangeables meme lorsque les moyennes sont egales, et la permutation
%   d'une statistique non studentisee derive. Permuter une statistique
%   studentisee de type Welch rend le test asymptotiquement valide sous
%   variances inegales tout en le gardant exact sous echangeabilite complete.
%   'permwelch' est donc la seule des trois variantes qui remplace a la fois le
%   test de normalite et le test d'homogeneite.
%   REFERENCES A VERIFIER AVANT CITATION : Janssen (1997), Studentized
%   permutation tests for non-iid hypotheses and the generalized Behrens-Fisher
%   problem ; Chung et Romano (2013), Exact and asymptotically robust
%   permutation tests. Ces deux entrees sont donnees de memoire et n'ont pas
%   ete verifiees sur le document original.
%
% La permutation ne change QUE la p-value. F standard, eta2, omega2 et donc
% l'ICC du verrou 2 sont calcules sur les donnees observees dans tous les cas.
if nargin<4 || isempty(diagMode),  diagMode='residual';  end
if nargin<5 || isempty(gate1Test), gate1Test='branching'; end
if nargin<6, PI=[]; end
[G,m]=size(Xp); N=G*m;
xg=mean(Xp,2); xx=mean(Xp(:));
MSb = m/(G-1)*sum((xg-xx).^2);
MSw = sum(sum((Xp-xg).^2))/(G*(m-1));
SSb=(G-1)*MSb; SSw=G*(m-1)*MSw; SSt=SSb+SSw;
eta2=SSb/SSt; omega2=(SSb-(G-1)*MSw)/(SSt+MSw);
Fstd = MSb/MSw;                          % F standard : ICC / effect sizes

% --- branche permutation, sans routage conditionnel ---------------------------
if ~strcmpi(gate1Test,'branching')
    if isempty(PI)
        error('anova_gate_paper:noPI', ...
              'La matrice de permutations PI est requise pour gate1Test = %s.', gate1Test);
    end
    switch lower(gate1Test)
        case 'perm',      unitSpec = true;  stat = 'F';     branch = 3;
        case 'permobs',   unitSpec = false; stat = 'F';     branch = 4;
        case 'permwelch', unitSpec = true;  stat = 'welch'; branch = 5;
        case 'permobswelch', unitSpec = false; stat = 'welch'; branch = 6;
        otherwise
            error('gate1Test inconnu : %s', gate1Test);
    end
    [p, Fw] = perm_pvalue(Xp, n, R, G, m, Fstd, PI, unitSpec, stat);
    if isempty(Fw), Fv = Fstd; else, Fv = Fw; end
    return
end

% --- cible des diagnostics d'hypothese selon le mode -------------------------
switch lower(diagMode)
    case 'residual'
        D = zeros(G,m);
        for g=1:G
            Xc = reshape(Xp(g,:), n, R);        % [specimens x ROIs] (col-major)
            D(g,:) = reshape(Xc - mean(Xc,2), 1, m);
        end
        dfw = G*n*(R-1);                        % ddl residuels (n moyennes/cond retirees)
    case 'pooled'
        D   = Xp;
        dfw = N-G;
    otherwise
        error('diagMode inconnu : %s (attendu ''residual'' ou ''pooled'')', diagMode);
end

normal_ok=true;
if exist('swtest','file')==2
    pv=zeros(G,1); for g=1:G,[~,pv(g)]=swtest(D(g,:)'); end
    normal_ok=all(pv>0.05/G);   % Bonferroni : sinon FWER = 1-0.95^G ~ 23% sous H0
end
Z=abs(D-median(D,2)); zg=mean(Z,2); zz=mean(Z(:));
Fbf=(m/(G-1)*sum((zg-zz).^2))/(sum(sum((Z-zg).^2))/dfw);
homo_ok = (1-fcdf(Fbf,G-1,dfw))>0.05;

Fstd = MSb/MSw;                          % rappel : calcule plus haut
branch = 0;                              % 0=F standard, 1=Welch, 2=Kruskal-Wallis
if ~normal_ok
    grp=repelem((1:G)',m); Xt=Xp.'; p=kruskalwallis(Xt(:),grp,'off'); Fv=Fstd;
    branch = 2;
elseif ~homo_ok
    s2=var(Xp,0,2); wg=m./s2; sw=sum(wg); xw=sum(wg.*xg)/sw;
    A=sum(wg.*(xg-xw).^2)/(G-1);
    L=sum((1-wg/sw).^2/(m-1));
    Fv=A/(1+2*(G-2)/(G^2-1)*L); nu=1/(3/(G^2-1)*L);
    p=1-fcdf(Fv,G-1,nu);
    branch = 1;
else
    Fv=Fstd; p=1-fcdf(Fstd,G-1,G*(m-1));
end
end
% -----------------------------------------------------------------------------
function [p, Fw_obs] = perm_pvalue(Xp, n, R, G, m, Fobs, PI, unitSpec, stat)
%PERM_PVALUE  p-value d'ANOVA par permutation, calcul vectorise.
%
%   UNITE DE PERMUTATION. Sous H0, ce sont les SPECIMENS qui sont echangeables
%   entre conditions : l'etiquette de condition s'attache au specimen, et les R
%   regions mesurees sur un meme specimen partagent son effet aleatoire. Elles
%   ne sont donc pas echangeables une a une. unitSpec = true permute les
%   specimens en bloc, avec leurs R regions ; unitSpec = false permute les
%   observations individuellement et n'est fourni que pour mesurer l'ecart.
%
%   STATISTIQUE. stat = 'F' permute le rapport des carres moyens, stat =
%   'welch' permute la statistique de Welch. La seconde reste valide sous
%   variances inegales, la premiere non (voir l'en-tete d'anova_gate_paper).
%
%   IDENTITE UTILISEE POUR stat = 'F'. SS_total est invariant par permutation,
%   donc SS_w = SS_t - SS_b et
%       F = [SS_b/(G-1)] / [(SS_t - SS_b)/(G(m-1))]
%   est strictement croissant en SS_b. Compter les F* superieurs ou egaux a
%   F_obs revient donc exactement a compter les SS_b* superieurs ou egaux a
%   SS_b,obs. Or SS_b ne depend que des sommes par condition. Aucune
%   approximation : le test est identique a celui obtenu en recalculant F a
%   chaque permutation.
%
%   Pour stat = 'welch', la statistique n'est pas monotone en SS_b et doit etre
%   recalculee. Elle ne demande cependant que la moyenne et la variance par
%   condition, donc les sommes ET les sommes de carres par condition. Ces deux
%   quantites sont additives sur les unites permutees, la meme vectorisation
%   s'applique et le cout reste negligeable.
%
%   CONVENTION. p = (1 + #{stat* >= stat_obs}) / (1 + B_perm) : l'unite ajoutee
%   compte l'etiquetage observe parmi les permutations possibles, ce qui
%   interdit une p-value nulle et garde le test exact (Phipson et Smyth, 2010).
%
%   PI est la matrice [B_perm x U] des permutations, U valant G*n au niveau
%   specimen et G*m au niveau observation. Elle est tiree une seule fois par
%   appel du pipeline et partagee par toutes les paires (s,k) : cela economise
%   le tirage et preserve la structure de dependance entre tests, ce qui est
%   coherent avec la correction FDR appliquee ensuite.
if nargin<9 || isempty(stat), stat='F'; end
B_perm = size(PI,1);
Fw_obs = [];

if unitSpec
    % sommes et sommes de carres par specimen, condition 1 puis 2, etc.
    Nspec = G*n;
    Sall  = zeros(Nspec, R);
    for g = 1:G
        Sall((g-1)*n + (1:n), :) = reshape(Xp(g,:), n, R);
    end
    s1 = sum(Sall, 2); s2 = sum(Sall.^2, 2);
    nu = n;                    % unites par condition
else
    % Xp est [G x m] : Xp(:) parcourt en COLONNES et entrelacerait les
    % conditions. La transposition est indispensable pour que les m premieres
    % valeurs appartiennent bien a la condition 1, etc.
    v  = reshape(Xp.', [], 1);
    s1 = v; s2 = v.^2;
    nu = m;
end
U = numel(s1);

% sommes par condition, affectation d'origine puis permutations
GS1o = sum(reshape(s1, nu, G), 1).';                       % [G x 1]
GS1  = squeeze(sum(reshape(s1(PI.'), nu, G, B_perm), 1));  % [G x B_perm]

switch lower(stat)
    case 'f'
        T = sum(s1);
        SSb_obs = sum(GS1o.^2)/m - T^2/(G*m);
        SSb     = sum(GS1.^2, 1)/m - T^2/(G*m);
        cnt = sum(SSb >= SSb_obs);
    case 'welch'
        GS2o = sum(reshape(s2, nu, G), 1).';
        GS2  = squeeze(sum(reshape(s2(PI.'), nu, G, B_perm), 1));
        Fw_obs = welch_stat(GS1o, GS2o, m, G);
        Fw     = welch_stat(GS1,  GS2,  m, G);
        cnt = sum(Fw >= Fw_obs);
    otherwise
        error('perm_pvalue:stat', 'stat inconnu : %s', stat);
end
p = (1 + cnt)/(1 + B_perm);
end
% -----------------------------------------------------------------------------
function Fw = welch_stat(GS1, GS2, m, G)
% Statistique de Welch a partir des sommes et sommes de carres par condition.
% GS1 et GS2 sont [G x B] : le calcul traite les B permutations d'un coup.
% Formule identique a celle de la branche Welch du logigramme.
xg = GS1/m;                               % moyennes par condition
vg = (GS2 - GS1.^2/m)/(m-1);              % variances non biaisees
w  = m./vg;
sw = sum(w, 1);
xw = sum(w.*xg, 1)./sw;
A  = sum(w.*(xg-xw).^2, 1)/(G-1);
L  = sum((1 - w./sw).^2/(m-1), 1);
Fw = A./(1 + 2*(G-2)/(G^2-1)*L);
end
% -----------------------------------------------------------------------------
function p_adj = bh_fdr(p)
M=numel(p); [ps,ix]=sort(p(:)); adj=ps.*M./(1:M)';
adj=flipud(cummin(flipud(adj)));
p_adj=zeros(M,1); p_adj(ix)=min(adj,1); p_adj=reshape(p_adj,size(p));
end
% -----------------------------------------------------------------------------
function [R2m_q, dAIC] = lme_quad(Xsk, Y)
%LME_QUAD  Ajuste Y ~ X + X^2 + (1|specimen) et compare au modele lineaire.
%   Renvoie le R2 marginal du modele quadratique et dAIC = AIC_lin - AIC_quad,
%   positif quand le terme quadratique merite son cout en parametres.
%
%   Les deux modeles sont ajustes par MAXIMUM DE VRAISEMBLANCE et non par REML :
%   ils different par leur structure d'effets FIXES, et les vraisemblances REML
%   de deux structures fixes differentes ne sont pas comparables.
[N,R] = size(Xsk);
x = Xsk(:); x = (x-mean(x))/std(x); y = Y(:);
sid = repelem((1:N)',R,1);
T = table(y, x, x.^2, nominal(sid), 'VariableNames',{'Y','X','X2','Sample'});
lin  = fitlme(T,'Y ~ X + (1|Sample)',      'FitMethod','ML');
quad = fitlme(T,'Y ~ X + X2 + (1|Sample)', 'FitMethod','ML');
b    = quad.Coefficients.Estimate;          % [b0 ; b1 ; b2]
varF = var(b(2)*x + b(3)*x.^2);             % part fixe, eq. R2 marginal
vc   = covarianceParameters(quad); varU = vc{1}(1,1);
R2m_q = varF/(varF + varU + quad.MSE);
dAIC  = lin.ModelCriterion.AIC - quad.ModelCriterion.AIC;
end
% -----------------------------------------------------------------------------
function [R2m,b1,se1] = lme_r2_marginal(Xsk, Y)
[N,R]=size(Xsk);
x=Xsk(:); x=(x-mean(x))/std(x); y=Y(:);
sid=repelem((1:N)',R,1);
T=table(y,x,nominal(sid),'VariableNames',{'Y','X','Sample'});
lme=fitlme(T,'Y ~ X + (1|Sample)','FitMethod','REML');
b1=lme.Coefficients.Estimate(2);
se1=lme.Coefficients.SE(2);                       % ecart-type a echantillon fixe
varF=var(b1*x);                                   % sigma2_f (eq 18)
vc=covarianceParameters(lme); varU=vc{1}(1,1);
varE=lme.MSE;
R2m=varF/(varF+varU+varE);
end
% -----------------------------------------------------------------------------
function [C, b1boot, info] = boot_C(Xsk, Y, b1_full, P, B_boot)
% Bootstrap specimen-level STRATIFIE par condition (n tirages dans chaque g).
% Renvoie la coherence de signe C, le VECTEUR des B_boot pentes, et le compte
% des ajustements aboutis.
%
% DENOMINATEUR. Un ajustement non abouti laisse son replicat a NaN et sort du
% denominateur. Compte au denominateur, il aurait ete lu comme une discordance
% de signe et aurait abaisse C sans raison.
%
% STANDARDISATION. Le centrage-reduction de X est FIXE sur l'echantillon
% complet. Recalcule a chaque replicat, il exprimerait chaque pente dans une
% unite differente et interdirait d'en tirer un ecart-type ou un intervalle.
% Le signe etant invariant par division par un scalaire positif, C est
% inchange par rapport a la version precedente.
[N,R] = size(Xsk); G = P.G; n = P.n;
mu_x = mean(Xsk(:)); sd_x = std(Xsk(:));
sgn = sign(b1_full);
b1boot = nan(B_boot,1);
for b = 1:B_boot
    idx = zeros(N,1);
    for g = 1:G
        base = (g-1)*n;
        idx(base+(1:n)) = base + randi(n,n,1);
    end
    try
        b1boot(b) = lme_b1_only(Xsk(idx,:), Y(idx,:), mu_x, sd_x);
    catch
        % ajustement non abouti : le replicat reste NaN
    end
end
ok   = ~isnan(b1boot);
info = struct('B_requested',B_boot,'B_valid',nnz(ok),'B_failed',nnz(~ok));
if info.B_valid == 0
    C = NaN;
else
    C = sum(sign(b1boot(ok)) == sgn) / info.B_valid;
end
end
function b1 = lme_b1_only(Xsk, Y, mu_x, sd_x)
% Ajustement REML, identique a celui de lme_r2_marginal : la distribution
% bootstrap doit porter sur le meme estimateur que l'estimation ponctuelle.
[N,R] = size(Xsk);
x = (Xsk(:)-mu_x)/sd_x; y = Y(:);
sid = repelem((1:N)',R,1);   % specimens reechantillonnes = clusters distincts
T = table(y,x,nominal(sid),'VariableNames',{'Y','X','Sample'});
lme = fitlme(T,'Y ~ X + (1|Sample)','FitMethod','REML');
b1 = lme.Coefficients.Estimate(2);
end
% -----------------------------------------------------------------------------
function ci = bca_ci(v, b1_full, Xsk, Y, alpha)
%BCA_CI  Intervalle bootstrap corrige du biais et accelere (Efron 1987).
%
%   POURQUOI PAS LE PERCENTILE. Sur ce cadre, la distribution bootstrap de la
%   pente est decalee vers zero pour tout descripteur qui n'est qu'un proxy de
%   la grandeur portant la relation : le reechantillonnage avec remise reduit
%   la variation intra-condition de X, et la pente s'attenue. Le decalage
%   atteint 0.86 ecart-type et il est reproductible sur replicats, ce qui
%   disqualifie l'intervalle percentile et justifie une correction.
%
%   DEUX CORRECTIONS.
%     z0     correction de biais, lue sur la fraction des pentes bootstrap
%            inferieures a l'estimation ponctuelle. Nulle si la distribution
%            est centree sur l'estimation.
%     accel  acceleration, qui corrige la variation de l'ecart-type le long de
%            l'echelle du parametre. Estimee par jackknife au niveau SPECIMEN,
%            c'est-a-dire a l'unite de reechantillonnage du bootstrap.
%
%   DOMAINE DE VALIDITE. La correction de biais deplace les quantiles demandes
%   vers les queues. Mesure faite sur 60 000 tirages : pour z0 autour de 0.28
%   la borne haute est convergee des B = 5000 ; pour z0 autour de 0.92 elle
%   derive encore a B = 20 000 ; pour z0 autour de 1.64 elle vaut le maximum de
%   l'echantillon quel que soit B, aucun tirage ne tombant au-dela du quantile
%   demande meme sur 60 000. BCa n'est donc applique que si les deux quantiles
%   corriges laissent au moins MIN_TAIL tirages de chaque cote. Sinon repli sur
%   le percentile, signale par out.b1_method.
n_b = numel(v);
MIN_TAIL = 10;

% --- correction de biais ------------------------------------------------------
% La proportion est bornee a 1/(2 n_b) de chaque cote : sans cela, une
% distribution entierement d'un cote de l'estimation renverrait z0 infini.
prop = sum(v < b1_full)/n_b;
prop = min(max(prop, 1/(2*n_b)), 1 - 1/(2*n_b));
z0   = norminv(prop);

% --- acceleration par jackknife specimen --------------------------------------
% Standardisation FIXE sur l'echantillon complet, identique a celle du
% bootstrap : les pentes jackknife doivent etre dans la meme unite.
[N,~] = size(Xsk);
mu_x = mean(Xsk(:)); sd_x = std(Xsk(:));
th = nan(N,1);
for i = 1:N
    keep = true(N,1); keep(i) = false;
    try
        th(i) = lme_b1_only(Xsk(keep,:), Y(keep,:), mu_x, sd_x);
    catch
        % ajustement non abouti : le specimen sort du calcul
    end
end
th = th(~isnan(th));
if numel(th) < 3
    accel = 0;
else
    d   = mean(th) - th;
    den = 6*(sum(d.^2))^1.5;
    if den < eps, accel = 0; else, accel = sum(d.^3)/den; end
end

% --- quantiles corriges -------------------------------------------------------
za = norminv(alpha/2); zb = norminv(1-alpha/2);
a1 = bca_alpha(z0, accel, za);
a2 = bca_alpha(z0, accel, zb);
qp = quantile(v, [alpha/2, 1-alpha/2]);
ci = struct('pct_lo',qp(1),'pct_hi',qp(2),'z0',z0,'accel',accel,'method','BCa');
if isnan(a1) || isnan(a2) || a1 >= a2 || a1*n_b < MIN_TAIL || (1-a2)*n_b < MIN_TAIL
    % Correction non definie, ou queues trop peu peuplees pour que les
    % quantiles corriges soient estimables a ce budget de tirages. Repli sur le
    % percentile brut. z0 et l'acceleration restent renvoyes : ils sont
    % informatifs en eux-memes, meme quand l'intervalle ne les utilise pas.
    ci.lo = qp(1); ci.hi = qp(2); ci.method = 'percentile';
else
    q = quantile(v, [a1 a2]);
    ci.lo = q(1); ci.hi = q(2);
end
end
function a = bca_alpha(z0, accel, z)
den = 1 - accel*(z0 + z);
if den <= 0
    a = NaN;
else
    a = normcdf(z0 + (z0 + z)/den);
end
end
% -----------------------------------------------------------------------------
function fragile = loso_flags(X, Y, out, P)
[N,R,S,K]=size(X); G=P.G; n=P.n; m0=(n-1)*R;
sampleCondition=repelem(1:G,n);
fragile=false(S,K);
for i=1:N
    keep=true(N,1); keep(i)=false;
    for s=1:S, for k=1:K
        y=[]; grp=[];
        for g=1:G
            sel=find(sampleCondition'==g & keep);
            v=reshape(X(sel,:,s,k),[],1);
            y=[y;v]; grp=[grp; g*ones(numel(v),1)]; %#ok<AGROW>
        end
        [p,tbl]=anova1(y,grp,'off');
        Fl=tbl{2,5}; e2=tbl{2,2}/tbl{4,2};
        iccl=max((Fl-1)/(Fl+m0-1),0);
        g2 = (p<P.q) & (e2>=P.eta2_thr) & (iccl>=P.icc_thr);
        if g2 ~= (out.decision(s,k)>=2), fragile(s,k)=true; end
    end, end
end
end
