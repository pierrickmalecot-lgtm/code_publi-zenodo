% version v10, 27/09/2026 : option P.noramp (defaut false) : amplitude de
%                          l'effet condition constante sur la bande de D1/D2,
%                          egale a sa valeur en s* ; test de la rampe
% version v9, 10/09/2026 : ajout du scenario 'latent' (aucun descripteur exact)
% version v8, 04/09/2026 : ajout du scenario 'ushape' (optimum interieur)
% version v7, 02/08/2026 : ajout des scenarios 'autocorr' (bruit ROI AR(1)),
% 'curvature' (reponse saturante) et 'randslope' (pente aleatoire)
% version v6, 31/07/2026 : panel de reference a cinq descripteurs, le correle mais peu fiable passe en scenario dedie
function [X, Y] = gen_volet1(P, scenario, rs)
%GEN_VOLET1  Generateur du volet 1 (rng externe), etendu aux scenarios.
%   VERSION 2 : ajout des scenarios de robustesse (violations d'hypotheses)
%   pour exercer les branches Welch / Kruskal-Wallis du pipeline.
%   VERSION 3 : diversification du panel de descripteurs.
%     D3 / D4 : bandes DISJOINTES (3-5 et 6-8) au lieu d'une bande commune
%               3-8. Teste la localisation de la bande et non la seule
%               detection ; supprime le recouvrement avec les deux bandes
%               fonctionnelles du scenario 'twoscale'.
%     D4      : bruit ROI seul (pas d'effet specimen) -> ICC attendu ~ 0.
%               C'est le null de reference pour le taux d'erreur de type I.
%     D5      : bruit specimen + ROI (structure hierarchique sans effet
%               condition). ATTENTION : l'effet specimen gonfle MSbetween
%               dans le pooling n*R, donc D5 peut presenter un F > 1 sous
%               l'hypothese nulle (pseudo-replication au sens de Hurlbert).
%               Ce n'est donc PAS un null de reference.
%   VERSION 6 : le panel de reference compte CINQ descripteurs, D1 a D5. Un
%   sixieme, correle a la reponse mais trop peu fiable pour fonder une
%   specification, n'existe que sous le scenario 'unreliable'. Il est construit
%   par un terme local partage avec la reponse, seul mecanisme qui separe
%   fiabilite et correlation : toute variabilite propre au descripteur fait
%   chuter le R2 plus vite que l'ICC, la pente s'attenuant et le residu
%   grossissant en meme temps.
%
%   scenario : 'nominal' | 'null' | 'twoscale' | 'unreliable'
%              | 'nonnormal' | 'heteroscedastic' | 'nonlinear'
%              | 'autocorr' | 'curvature' | 'randslope' | 'ushape'
%              | 'latent'
%   rs : RandStream (ou [] pour rng global, reproduit rng(42) du volet 1)
%   X : [N x R x S x K] ; Y : [N x R]
%
%   Scenarios v2 :
%     'nonnormal'       bruit intra log-normal (asymetrique) -> Shapiro-Wilk
%                       rejette par condition -> branche Kruskal-Wallis.
%     'heteroscedastic' sigma_within croissant avec la condition -> Brown-
%                       Forsythe rejette -> branche Welch (~84% Welch, ~16% KW
%                       sur 50 tirages ; les deux sont robustes, l'objectif est
%                       de sortir de la branche F standard).
%   Scenarios v7 (section 6.4) :
%     'autocorr'        bruit intra AR(1) le long des ROI, coefficient
%                       P.rho_spat -> dependance entre regions de mesure, la
%                       variance marginale restant sigma_within^2
%     'curvature'       reponse saturante de taux P.curv -> relation monotone
%                       mais courbe, le R2 lineaire se degrade quand le
%                       Spearman ne bouge pas
%
%     'nonlinear'       Y = fonction monotone saturante de X -> R2 lineaire
%                       sous-estime l'association que le Spearman revele.
%
if isempty(rs), rfun = @(varargin) randn(varargin{:});
else,           rfun = @(varargin) randn(rs, varargin{:}); end
G=P.G; n=P.n; R=P.R; S=P.S; K=P.K; N=P.N;
sampleCondition = repelem(1:G, n);
X = zeros(N,R,S,K);
Z = zeros(N,R);                 % terme local, scenario 'unreliable' seulement
% v6 : le panel de reference compte CINQ descripteurs. Le descripteur correle
% mais peu fiable a ete sorti du panel et fait l'objet d'un scenario dedie,
% 'unreliable', ou il occupe la sixieme colonne. Le terme local qui le
% construit n'a pas de contrepartie physique mesurable sur une surface reelle :
% il sert a fabriquer une situation, au meme titre que le tirage independant de
% la reponse dans le scenario 'null'. Le laisser dans le panel de reference
% aurait abaisse tous les R2 du jeu principal sans raison defendable.
isUnrel = strcmp(scenario,'unreliable');

% bandes selon scenario (D1/D2 pilotes)
if strcmp(scenario,'twoscale'), b12 = {[2 4],[7 9]}; else, b12 = {[5 7],[5 7]}; end
% v5 : un seul discriminant seul, D3, sur la bande basse
b3 = [3 5];

% --- v2 : parametres de violation d'hypothese (bruit intra w) -----------------
noiseMode = 'normal';
if strcmp(scenario,'nonnormal'),       noiseMode = 'lognormal'; end
if strcmp(scenario,'heteroscedastic'), noiseMode = 'hetero';    end
if strcmp(scenario,'autocorr'),        noiseMode = 'ar1';       end
sg_ln = 0.6;                                   % log-normale
mln   = exp(sg_ln^2/2);                         % E[exp(sg Z)]
vln   = (exp(sg_ln^2)-1)*exp(sg_ln^2);          % Var[exp(sg Z)]

for i=1:N
    g = sampleCondition(i); mu_g = P.condLevels(g);
    a_i = P.sigma_between * rfun(1);
    % v3 : effet specimen INDEPENDANT pour D6. Le a_i commun est partage par
    % D1 (qui pilote Y), donc l'utiliser pour D6 creerait une correlation
    % reelle entre D6 et Y via l'effet specimen. D6 doit rester un null
    % propre vis-a-vis de la reponse tout en portant une structure
    % hierarchique.
    a_i6 = P.sigma_between * rfun(1);
    % v7 : etat AR(1) le long des ROI, scenario 'autocorr'. L'ordre des
    % boucles etant i > r > s > k, la dependance le long de r doit etre
    % portee par un tampon indexe (s,k), remis a zero a chaque specimen :
    % les ROI d'un specimen sont correlees entre elles, pas d'un specimen
    % a l'autre.
    wPrev = zeros(S,K);
    for r=1:R
        if isUnrel, z_ir = P.sigma_z * rfun(1); else, z_ir = 0; end
        Z(i,r) = z_ir;
        for s=1:S
            for k=1:K
                % --- bruit intra w selon le scenario (v2) ---
                switch noiseMode
                    case 'lognormal'
                        z = rfun(1);
                        w = P.sigma_within * (exp(sg_ln*z) - mln) / sqrt(vln);
                    case 'hetero'
                        hetScale = 0.4 + 1.2*(g-1)/(G-1);   % ratio variance ~16x
                        w = P.sigma_within * hetScale * rfun(1);
                    case 'ar1'
                        % AR(1) stationnaire le long des ROI. La variance
                        % marginale reste sigma_within^2 : seule la dependance
                        % entre ROI voisines change, les composantes de
                        % variance restent comparables au cas nominal, et
                        % l'effet observe est donc imputable a la seule
                        % violation d'independance.
                        if r==1
                            w = P.sigma_within * rfun(1);
                        else
                            w = P.rho_spat*wPrev(s,k) ...
                                + P.sigma_within*sqrt(1-P.rho_spat^2)*rfun(1);
                        end
                        wPrev(s,k) = w;
                    otherwise
                        w = P.sigma_within * rfun(1);
                end
                switch k
                    case {1,2}                      % pilotes (rampe volet 1)
                        bd = b12{k};
                        if s>=bd(1) && s<=bd(2)
                            % Rampe de reference : amplitude x1.0, x1.3, x1.6
                            % sur les echelles 5, 6, 7. Avec P.noramp = true,
                            % amplitude constante x1.3, valeur en s* = 6 : D1
                            % en s* et donc Y sont inchanges.
                            if isfield(P,'noramp') && P.noramp
                                signal = mu_g * 1.3;
                            else
                                signal = mu_g * (1 + 0.3*(s-bd(1)));
                            end
                        else
                            signal = 0.5 + 0.05*rfun(1);
                        end
                        X(i,r,s,k) = signal + a_i + w;
                    case 3                          % discriminant seul (v5)
                        % Signal condition sur la bande basse 3 a 5, absent de
                        % l'equation de la reponse. Sa bande differe de celle
                        % de D1 et D2, ce qui verifie que la procedure attribue
                        % chaque bande au bon descripteur et ne se contente pas
                        % de detecter qu'il y a du signal quelque part.
                        if s>=b3(1) && s<=b3(2)
                            signal = mu_g*0.8 + 0.2*sin(2*pi*mu_g*3);
                        else
                            signal = 0.5 + 0.05*rfun(1);
                        end
                        X(i,r,s,k) = signal + a_i + w;
                    case 4                          % bruit ROI SEUL (v3)
                        % aucun effet condition, aucun effet specimen :
                        % variance intra pure -> ICC attendu ~ 0.
                        X(i,r,s,k) = 0.5 + w;
                    case 5                          % bruit SPECIMEN + ROI (v3)
                        % aucun effet condition, mais structure hierarchique
                        % presente (a_i) : exerce le verrou 2 sur un cas
                        % nouveau (variance structuree mais non discriminante).
                        X(i,r,s,k) = 0.5 + a_i6 + w;
                    case 6                          % correle mais peu fiable (v6)
                        % N'EXISTE QUE sous le scenario 'unreliable'. Hors de
                        % ce scenario le terme local z est nul et cette colonne
                        % reproduit simplement D1, sans effet sur les analyses
                        % qui ne portent que sur les cinq premieres.
                        % Meme bande que D1, 5 a 7, sans rampe, plus un terme
                        % z_ir propre a la region de mesure et PARTAGE avec la
                        % reponse, ou il entre avec le coefficient P.gamma_z.
                        % Ce partage est le seul mecanisme qui separe les deux
                        % verrous. Toute variabilite propre a D4, quel que soit
                        % le niveau ou on la place, fait chuter le R2 plus vite
                        % que l'ICC : la pente s'attenue et le residu grossit,
                        % le R2 est donc frappe deux fois. Un terme partage, au
                        % contraire, gonfle la dispersion intra-condition, donc
                        % abaisse l'ICC, tout en etant predictif de la reponse,
                        % donc en soutenant le R2.
                        % Lecture physique : une heterogeneite locale de matiere
                        % qui affecte a la fois la topographie relevee a cet
                        % endroit et le comportement fonctionnel local. Ce
                        % descripteur est le seul que le verrou de fiabilite
                        % doive intercepter alors que le verrou de correlation
                        % l'aurait accepte, d'ou son scenario dedie.
                        if s>=5 && s<=7
                            signal = mu_g;
                        else
                            signal = 0.5 + 0.05*rfun(1);
                        end
                        X(i,r,s,k) = signal + a_i + w + z_ir;
                end
            end
        end
    end
end

% ---- Y ----------------------------------------------------------------------
Y = zeros(N,R);
switch scenario
    case {'nominal','nonnormal','heteroscedastic','unreliable','autocorr'}
        % Y lineaire en D1 (les scenarios de robustesse ne changent que le
        % bruit intra de X, pas la relation fonctionnelle Y(X))
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                Y(i,r) = P.beta0 + P.beta1*X(i,r,P.s_star,1) + u_i + P.sigma_func*rfun(1) ...
                         + P.gamma_z * Z(i,r);   % nul hors scenario 'unreliable'
            end
        end
    case 'latent'
        % v9 : la reponse depend d'une variable LATENTE, dont aucun descripteur
        % n'est la mesure exacte.
        %
        % POURQUOI. Dans tous les autres scenarios, Y est construit directement
        % sur X(:,:,s_star,1) : il existe donc une paire pour laquelle la
        % relation est algebriquement exacte, sans terme d'erreur sur le
        % predicteur. Le reechantillonnage bootstrap laisse alors sa pente
        % inchangee et la correction de biais z0 s'y annule, ce qui produit le
        % contraste net observe entre D1 et D2 a l'echelle generatrice.
        %
        % Une surface reelle ne se comporte pas ainsi : aucun descripteur n'est
        % la cause du comportement fonctionnel, tous sont des indicateurs
        % partiels de la microgeometrie qui, elle, le produit.
        %
        % SENS DE CONSTRUCTION. Il est decisif et une premiere version s'y est
        % trompee. Poser T = X - e puis Y = f(T) revient a ecrire
        % Y = f(X) - f(e) : le bruit se retrouve dans la REPONSE, ce qui elargit
        % le residu sans attenuer la pente, et le descripteur reste le
        % predicteur exact. Il faut l'inverse : la variable latente est le
        % descripteur tel qu'il a ete engendre, et c'est le descripteur OBSERVE
        % qui en est une mesure bruitee.
        %     T          = X engendre, non observable
        %     Y          = beta0 + beta1*T + u_i + eps
        %     X observe := T + bruit de mesure
        % Le bruit de mesure est ajoute a TOUS les descripteurs et a toutes les
        % echelles, une erreur de mesure n'ayant pas de raison d'epargner les
        % autres paires. Cela affecte aussi les deux premieres portes, ce qui
        % est voulu : un descripteur mesure avec erreur discrimine moins bien.
        %
        % P.lat_noise fixe l'ecart-type de ce bruit, en fraction de l'ecart-type
        % du descripteur porteur. A 0 le scenario redevient identique a
        % 'nominal'.
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                Y(i,r) = P.beta0 + P.beta1*X(i,r,P.s_star,1) ...
                         + u_i + P.sigma_func*rfun(1);
            end
        end
        % La reponse etant fixee sur les valeurs engendrees, le bruit de mesure
        % est ajoute ensuite : aucune paire n'est plus le predicteur exact.
        if P.lat_noise > 0
            sdX = std(reshape(X(:,:,P.s_star,1),[],1));
            for i=1:N
                for r=1:R
                    for s=1:S
                        for k=1:K
                            X(i,r,s,k) = X(i,r,s,k) + P.lat_noise*sdX*rfun(1);
                        end
                    end
                end
            end
        end
    case 'randslope'
        % v7 : pente aleatoire. Y = beta0 + (beta1 + u1_i) X + u0_i + eps.
        % Sert a verifier que la regle de selection de modele (AIC) retient la
        % specification complete quand sigma_u1 est non negligeable devant
        % beta1, et la specification reduite sinon. Avec P.sigma_u1 = 0 le
        % scenario est identique a 'nominal' quant a la structure de Y.
        for i=1:N
            u_i  = P.sigma_u  * rfun(1);
            u1_i = P.sigma_u1 * rfun(1);
            for r=1:R
                Y(i,r) = P.beta0 + (P.beta1 + u1_i)*X(i,r,P.s_star,1) ...
                         + u_i + P.sigma_func*rfun(1);
            end
        end
    case 'curvature'
        % v7 : reponse courbe mais monotone, taux de courbure P.curv.
        %
        % Un terme quadratique convexe ne convient pas : sur la plage positive
        % de X occupee par le plan, corr(x, x^2) vaut 0.97, donc x^2 est
        % quasi lineaire en x. L'ajouter ne cree aucune courbure exploitable
        % et augmente meme le R2 en raidissant la pente.
        %
        % La forme saturante, elle, brise la linearite tout en restant
        % strictement croissante : le R2 marginal se degrade avec P.curv alors
        % que le Spearman reste a 1. C'est la separation des deux indices
        % recherchee en section 6.2. La normalisation par (1-exp(-c)) garde la
        % meme amplitude de reponse quel que soit le taux, de sorte que seule
        % la forme change d'un niveau a l'autre.
        %
        % Repere sur le plan de reference : c = 1 -> R2 ~ 0.98, c = 3 -> 0.86,
        % c = 5 -> 0.72, Spearman = 1.00 dans les trois cas.
        cc = P.curv;
        Xs = X(:,:,P.s_star,1);
        xmin = min(Xs(:)); xrng = max(Xs(:)) - xmin;
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                xn = (X(i,r,P.s_star,1) - xmin) / xrng;
                if cc < 1e-6
                    g = xn;                                   % limite lineaire
                else
                    g = (1 - exp(-cc*xn)) / (1 - exp(-cc));   % saturante
                end
                Y(i,r) = P.beta0 + P.beta1*g + u_i + P.sigma_func*rfun(1);
            end
        end
    case 'ushape'
        % v8 : reponse a OPTIMUM INTERIEUR, en cloche. C'est le cas que la
        % porte 3 rejette en silence : la relation existe et est forte, mais
        % elle n'est ni lineaire ni monotone, donc le R2 marginal ET le
        % Spearman tombent ensemble.
        %
        % Forme : g(xn) = 1 - ((xn - x0)/s)^2 sur xn normalise dans [0,1],
        % avec s = max(x0, 1-x0) pour que g couvre [0,1] quelle que soit la
        % position de l'optimum. A x0 = 0.5 la parabole est symetrique sur la
        % plage et corr(xn, g) est nulle par symetrie : la pente lineaire
        % s'annule et le Spearman aussi. En deplacant x0 vers un bord, la
        % branche descendante domine, une composante monotone reapparait et le
        % Spearman remonte en valeur absolue.
        %
        % Contrairement au scenario 'curvature', qui reste strictement
        % croissant et ne degrade que le R2, celui-ci degrade les deux indices
        % a la fois. C'est cette conjonction, sur une paire ayant franchi les
        % deux premieres portes, qui constitue le signalement.
        x0 = P.ushape_x0;
        ss = max(x0, 1-x0);
        Xs = X(:,:,P.s_star,1);
        xmin = min(Xs(:)); xrng = max(Xs(:)) - xmin;
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                xn = (X(i,r,P.s_star,1) - xmin) / xrng;
                g  = 1 - ((xn - x0)/ss)^2;
                Y(i,r) = P.beta0 + P.beta1*g + u_i + P.sigma_func*rfun(1);
            end
        end
    case 'nonlinear'
        % relation monotone saturante : R2 lineaire sous-estime l'association
        % que le Spearman revele (role propre de Spearman, section 3.8.2)
        Xs = X(:,:,P.s_star,1);
        xmin = min(Xs(:)); xrng = max(Xs(:)) - xmin;
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                xn = (X(i,r,P.s_star,1) - xmin) / xrng;      % normalise [0,1]
                Y(i,r) = P.beta0 + P.beta1*(1 - exp(-5*xn)) + u_i + P.sigma_func*rfun(1);
            end
        end
    case 'null'
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                Y(i,r) = P.beta0 + u_i + P.sigma_func*rfun(1);
            end
        end
    case 'twoscale'
        XA = X(:,:,3,1); XB = X(:,:,8,2);
        T  = P.sigma_u^2 / 0.10;                 % parts 35/35/10 (voir 4.8)
        b1 = sqrt(0.35*T/var(XA(:))); b2 = sqrt(0.35*T/var(XB(:)));
        for i=1:N
            u_i = P.sigma_u * rfun(1);
            for r=1:R
                Y(i,r) = P.beta0 + b1*XA(i,r) + b2*XB(i,r) + u_i + P.sigma_func*rfun(1);
            end
        end
end
end
