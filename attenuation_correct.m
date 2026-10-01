% version v3, 17/09/2026 : intervalle bootstrap recalculant kappa a chaque tirage,
%                          domaine de validite mis a jour au seuil 0.50
% version v2, 11/09/2026 : notation kappa (lambda etant pris), seuil porte a 0.50
% version v1, 05/09/2026 : correction d'attenuation pour la configuration (b)
function out = attenuation_correct(Xsk, Y, alpha, grp, B, rs)
%ATTENUATION_CORRECT  Pente corrigee de l'erreur sur le descripteur moyen,
%   configuration (b) : une seule mesure fonctionnelle par specimen.
%
%   LE PROBLEME. En configuration (b), la reponse Y_i n'existe qu'au niveau du
%   specimen et la correlation se fait par moindres carres sur le descripteur
%   moyenne sur les R regions. Or cette moyenne est elle-meme entachee d'erreur,
%   R regions ne caracterisant pas entierement une surface. Un predicteur
%   entache d'erreur tire la pente estimee vers zero. La section sur les
%   variantes de reponse signale cette attenuation et recommande de rapporter la
%   dispersion entre regions ; cette fonction la quantifie et la corrige.
%
%   CE QU'IL FAUT ESTIMER. La fiabilite du predicteur, soit la part de la
%   variance de X_barre qui est du vrai signal specimen et non de l'erreur de
%   moyennage :
%       kappa = (var(X_barre) - MSW/R) / var(X_barre)
%   ou MSW est le carre moyen intra-specimen de X, donc MSW/R la variance
%   d'erreur de la moyenne de R regions. La pente non attenuee s'obtient en
%   divisant la pente observee par kappa, et le R2 de meme.
%
%   NOTATION. kappa et non lambda : lambda designe deja, section 3.9.3, le
%   rapport de l'ecart-type bootstrap a l'ecart-type du modele. Les deux sont
%   des rapports lies a l'incertitude, ce qui rend la distinction d'autant plus
%   necessaire.
%
%   ATTENTION A L'ECHELLE. La variance a mettre au denominateur est celle de
%   X_barre SUR L'ECHANTILLON DE REGRESSION, conditions confondues, puisque la
%   regression porte sur les N = G*n specimens. Une fiabilite calculee a
%   l'interieur des conditions donnerait une valeur bien plus basse et
%   surcorrigerait lourdement. Ce n'est pas non plus l'ICC de la porte 2, qui
%   est construit sur le F inter-conditions et repond a une autre question.
%
%   DOMAINE DE VALIDITE, MESURE. Plan de reference, G=5, n=5, R=5, 999
%   replicats par point, KAPPA_MIN = 0.50 (run_attenuation, v3).
%     attenuation      pente non corrigee   corrigee, mediane   corrigee, moyenne
%     2 a 33 %         -2 a -32 %           -0.5 a +0.5 %       0 a +2.5 %
%     41 %             -40 %                -7 %                -5 %
%   Jusqu'a une attenuation d'un tiers, la correction retablit la pente. Pres
%   du seuil (kappa vrai 0.59), elle se degrade : 26 pour cent des jeux ne sont
%   plus corriges, et les autres sous-estiment la pente d'environ 7 pour cent,
%   le seuil ecartant de preference les tirages ou kappa estime, donc la
%   correction, aurait ete le plus fort.
%   Configuration (c), T essais par specimen (run_attenuation_c, v1) : memes
%   ordres de grandeur, avec T = 3 pour tous comme avec T tire dans {2, 3, 4}.
%   A T constant, la pente du modele mixte est exactement celle des moindres
%   carres sur les moyennes de specimen.
%
%   KAPPA_MIN est fixe a 0.50. A 0.30, valeur essayee d'abord, la correction
%   s'emballait : la moyenne de la pente corrigee derivait de +7.6 pour cent a
%   41 pour cent d'attenuation, le quantile a 95 pour cent atteignait 6.25 pour
%   une pente vraie de 3.5, et le R2 corrige saturait a 1.000 : les kappa
%   estimes entre 0.30 et 0.50 etaient encore corriges. Le kappa estime peut
%   aussi devenir negatif, la variance d'erreur estimee depassant alors la
%   variance totale observee par simple fluctuation d'echantillonnage ; ces cas
%   sont ecartes quel que soit le seuil. Sous le seuil, la correction n'est pas
%   appliquee et la dispersion entre regions doit etre rapportee telle quelle.
%
%   ARGUMENTS.
%     Xsk    descripteur, [N x R] pour la paire (s,k) consideree.
%     Y      reponse fonctionnelle, un scalaire par specimen, [N x 1].
%     alpha  niveau pour l'intervalle. Defaut 0.05.
%     grp    condition de chaque specimen, [N x 1]. Facultatif.
%     B      nombre de tirages bootstrap. Defaut 0, pas de bootstrap.
%     rs     flux aleatoire pour le bootstrap. Facultatif.
%
%   INTERVALLE BOOTSTRAP. L'intervalle ci_corr divise l'ecart-type de la pente
%   par kappa estime, comme si kappa etait connu. Mesure par
%   run_attenuation_c, il ne couvre la pente vraie que dans 85 a 89 pour cent
%   des cas des que l'attenuation depasse 20 pour cent, pour 95 vises. Avec
%   B > 0, les specimens sont reechantillonnes avec remise, a l'interieur de
%   chaque condition si grp est fourni, et kappa est recalcule sur chaque
%   tirage. Sa variance entre dans l'ecart-type de la pente corrigee par
%   propagation au premier ordre :
%       se_boot^2 = se^2/kappa^2 + b1^2 * var_boot(kappa) / kappa^4
%   ci_boot = b1_corr +/- t * se_boot.
%   La variance bootstrap est multipliee par n/(n-1), n etant le nombre de
%   specimens par condition : a n = 5, le bootstrap stratifie sous-estime la
%   dispersion d'environ un cinquieme.
%   CHOIX MESURE (simulation Python, 800 jeux par point, plan de reference).
%   Couverture a 95 pour cent visee, de 2 a 41 pour cent d'attenuation :
%     intervalle ci_corr                         85 a 96 pour cent
%     percentiles de la pente corrigee bootstrap 86 a 97 pour cent, trop
%                                                etroit a faible attenuation
%     propagation, avec facteur n/(n-1)          95 a 97 pour cent
%     propagation, sans ce facteur               92 a 96 pour cent
%   Aucun seuil n'est applique aux tirages : un kappa de tirage proche de zero
%   elargit l'intervalle, ce qui est l'information recherchee. La part de ces
%   tirages est rapportee.
%
%   SORTIES, dans la structure out.
%     b1, b1_corr        pente observee et corrigee
%     R2, R2_corr        coefficient de determination, idem
%     kappa             fiabilite estimee du descripteur moyen
%     attenuation        1 - kappa, part de pente perdue
%     se, se_corr        ecart-type de la pente, avant et apres correction
%     ci_corr            intervalle de la pente corrigee
%     reliable           false si kappa est trop bas pour corriger
%     se_boot, ci_boot   ecart-type et intervalle de la pente corrigee,
%                        incertitude sur kappa comprise (B > 0)
%     var_kappa_boot     variance bootstrap de kappa, facteur n/(n-1) inclus
%     boot_low_kappa     part des tirages ou kappa est sous KAPPA_MIN
%     msw, var_xbar      quantites intermediaires, pour report
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4, grp = []; end
if nargin < 5 || isempty(B), B = 0; end
if nargin < 6 || isempty(rs), rs = RandStream.getGlobalStream; end
KAPPA_MIN = 0.50;   % en deca, la correction n'est plus fiable

[N, R] = size(Xsk);
Y = Y(:);
if numel(Y) ~= N
    error('attenuation_correct:size', ...
          'Y doit contenir une valeur par specimen (%d attendues, %d fournies).', ...
          N, numel(Y));
end

xbar = mean(Xsk, 2);

% --- regression observee ------------------------------------------------------
D  = [ones(N,1) xbar];
b  = D\Y;
r  = Y - D*b;
s2 = sum(r.^2)/(N-2);
C  = s2 * inv(D.'*D);                                     %#ok<MINV>
b1 = b(2);
se = sqrt(C(2,2));
R2 = 1 - sum(r.^2)/sum((Y-mean(Y)).^2);

% --- fiabilite du descripteur moyen -------------------------------------------
msw      = sum(sum((Xsk - xbar).^2)) / (N*(R-1));   % variance intra-specimen
var_err  = msw / R;                                  % erreur de la moyenne
var_xbar = var(xbar, 0);
kappa   = (var_xbar - var_err) / var_xbar;

out = struct('b1',b1,'se',se,'R2',R2,'kappa',kappa, ...
             'attenuation',1-kappa,'msw',msw,'var_xbar',var_xbar, ...
             'reliable',kappa >= KAPPA_MIN);

if ~out.reliable
    % Correction non appliquee : sous ce seuil le rapport est trop instable
    % pour etre interprete. La dispersion entre regions doit alors etre
    % rapportee telle quelle, et le plan revu.
    out.b1_corr = NaN; out.R2_corr = NaN; out.se_corr = NaN;
    out.ci_corr = [NaN NaN];
    out.se_boot = NaN; out.ci_boot = [NaN NaN]; out.boot_low_kappa = NaN;
    out.var_kappa_boot = NaN;
    return
end

% --- correction ---------------------------------------------------------------
% La division par kappa retire l'attenuation et elargit d'autant l'incertitude :
% l'intervalle corrige est plus large que l'observe, ce qui est le prix de la
% correction et doit etre rapporte avec elle.
out.b1_corr = b1 / kappa;
out.se_corr = se / kappa;
out.R2_corr = min(R2 / kappa, 1);
tq = tinv(1-alpha/2, N-2);
out.ci_corr = out.b1_corr + [-1 1]*tq*out.se_corr;

% --- bootstrap : kappa recalcule sur chaque tirage ---------------------------
out.se_boot = NaN; out.ci_boot = [NaN NaN]; out.boot_low_kappa = NaN;
out.var_kappa_boot = NaN;
if B > 0
    if isempty(grp), grp = ones(N,1); end
    grp = grp(:);
    ug  = unique(grp).';
    idx = zeros(B, N);
    pos = 0;
    for g = ug
        m  = find(grp == g);
        ng = numel(m);
        idx(:, pos+1:pos+ng) = m(randi(rs, ng, B, ng));
        pos = pos + ng;
    end
    ssw = sum((Xsk - xbar).^2, 2);                 % somme intra, par specimen
    xb  = xbar(idx);                               % [B x N]
    vxb = var(xb, 0, 2);
    kb  = (vxb - sum(ssw(idx), 2) / (N*(R-1)) / R) ./ vxb;
    nc  = N / numel(ug);                           % specimens par condition
    vk  = var(kb) * nc / (nc - 1);
    out.var_kappa_boot = vk;
    out.se_boot = sqrt(se^2 / kappa^2 + b1^2 * vk / kappa^4);
    out.ci_boot = out.b1_corr + [-1 1]*tq*out.se_boot;
    out.boot_low_kappa = mean(kb < KAPPA_MIN);
end
end
