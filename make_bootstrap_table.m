% version v1, 02/09/2026 : tableau de sortie du bootstrap, paires retenues
%MAKE_BOOTSTRAP_TABLE  Ce que le bootstrap apporte, une fois C mis de cote.
%
%   Usage : lancer run_volet1, puis
%              make_bootstrap_table(out)
%           ou make_bootstrap_table(out, 'bootstrap_table.csv') pour ecrire.
%
%   POURQUOI PAS C. La coherence de signe sature des que l'intervalle ne
%   traverse pas zero, ce qui est le cas de toute paire ayant franchi les deux
%   verrous amont. Sur le tirage de reference elle vaut 1.000 sur les neuf
%   paires retenues, alors que la precision relative varie d'un facteur quatre
%   entre ces memes paires. C est constante sur son domaine d'emploi.
%
%   TROIS QUANTITES, TROIS QUESTIONS DISTINCTES.
%
%   rp = ecart-type bootstrap / |pente|
%       Precision relative de la pente. Repond a : la relation est-elle estimee
%       assez finement pour servir de base a une specification. Sans dimension,
%       comparable entre descripteurs d'unites differentes.
%
%   ratio = ecart-type bootstrap / ecart-type parametrique
%       L'ecart-type rendu par fitlme est calcule a echantillon de specimens
%       FIXE : il ne compte que le bruit de mesure et la dispersion des
%       interceptes. Le bootstrap reechantillonne les specimens, il ajoute donc
%       la variabilite due au fait que ces n specimens ne sont qu'un tirage
%       parmi les possibles. Un ratio de 2 signifie que la moitie de
%       l'incertitude vient du choix des specimens : ajouter des regions de
%       mesure ne servirait a rien, il faut plus de specimens. Un ratio proche
%       de 1 signifie l'inverse.
%
%   z0 = correction de biais, norminv de la fraction des pentes bootstrap
%        inferieures a l'estimation
%       Le reechantillonnage avec remise reduit la variation intra-condition de
%       X. Pour le descripteur qui porte reellement la relation, cela ne coute
%       rien. Pour un descripteur qui n'en est qu'un proxy, cette variation est
%       precisement ce qui porte l'information propre a la paire, et la pente
%       s'attenue. z0 mesure cette attenuation, sans dimension. Sur le jeu de
%       reference il vaut 0.07 sur la paire porteuse et monte a 1.64 sur les
%       paires les plus eloignees d'elle.
%
%   METHODE. La colonne methode indique si l'intervalle est un BCa ou un repli
%   sur le percentile : au-dela de z0 voisin de 0.9, les quantiles corriges
%   tombent dans une queue trop peu peuplee pour etre estimables a B = 1000.
function T = make_bootstrap_table(out, csvfile)
if nargin < 2, csvfile = ''; end
if ~isfield(out,'b1_se')
    error('make_bootstrap_table:noBootstrap', ...
          'Lancer le pipeline avec opts.bootstrap = true.');
end

[si,ki] = find(out.decision==3);
if isempty(si), fprintf('Aucune paire retenue.\n'); T = table(); return; end
[~,ord] = sortrows([ki si]); si = si(ord); ki = ki(ord);

n = numel(si);
scale = si(:); desc = ki(:);
beta1 = nan(n,1); lo = beta1; hi = beta1; se_boot = beta1; se_par = beta1;
ratio = beta1; rp = beta1; z0 = beta1; C = beta1;
method = cell(n,1);
for j = 1:n
    s = si(j); k = ki(j);
    beta1(j)   = out.beta1(s,k);
    lo(j)      = out.b1_lo(s,k);
    hi(j)      = out.b1_hi(s,k);
    se_boot(j) = out.b1_se(s,k);
    se_par(j)  = out.b1_se_par(s,k);
    ratio(j)   = out.b1_ratio(s,k);
    rp(j)      = out.b1_rp(s,k);
    z0(j)      = out.b1_z0(s,k);
    C(j)       = out.C(s,k);
    method{j}  = out.b1_method{s,k};
end

T = table(desc, scale, beta1, lo, hi, se_boot, se_par, ratio, rp, z0, C, ...
    method, 'VariableNames', {'descriptor','scale','beta1','ci_lo','ci_hi', ...
    'se_boot','se_par','ratio','rp','z0','C','method'});

fprintf('\n=== Sorties du bootstrap, paires retenues ===\n');
fprintf('rp    = se_boot / |beta1|      precision relative de la pente\n');
fprintf('ratio = se_boot / se_par       part de l incertitude due aux specimens\n');
fprintf('z0    = correction de biais    attenuation propre aux descripteurs proxy\n\n');
fprintf(' desc  ech    beta1      intervalle 95%%     se_boot  se_par  ratio     rp     z0      C   methode\n');
for j = 1:n
    fprintf(' D%-3d %4d  %7.3f  [%6.3f %6.3f]  %7.3f %7.3f %6.2f %6.3f %6.2f %6.3f   %s\n', ...
        desc(j), scale(j), beta1(j), lo(j), hi(j), se_boot(j), se_par(j), ...
        ratio(j), rp(j), z0(j), C(j), method{j});
end

fprintf('\nEtendue de C     : %.3f a %.3f  (facteur %.2f)\n', ...
    min(C), max(C), max(C)/min(C));
fprintf('Etendue de rp    : %.3f a %.3f  (facteur %.2f)\n', ...
    min(rp), max(rp), max(rp)/min(rp));
fprintf('Etendue du ratio : %.2f a %.2f  (facteur %.2f)\n', ...
    min(ratio), max(ratio), max(ratio)/min(ratio));
fprintf('Etendue de z0    : %.2f a %.2f\n', min(z0), max(z0));

if ~isempty(csvfile)
    writetable(T, csvfile);
    fprintf('\nTableau ecrit dans %s\n', csvfile);
end
end
