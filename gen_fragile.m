% version v1, 31/07/2026 : injection d une anomalie sur un specimen
function [X, Y, meta] = gen_fragile(P, rs, opts)
%GEN_FRAGILE  Scenario "specimen aberrant" pour exercer le DRAPEAU FRAGILE du LOSO.
%   Part du dataset nominal (gen_volet1) et injecte une anomalie sur UN SEUL
%   specimen, a UNE SEULE paire (s,k) cible, de facon a rendre la decision de
%   verrou dependante de la presence de ce specimen : sur donnees completes la
%   paire echoue (ou passe) un verrou ; en retirant ce specimen, le statut
%   s'inverse -> loso_flags(...) leve fragile(s,k)=true.
%
%   Deux mecanismes physiquement interpretables (usinage) :
%     'variance'  un specimen aux ROIs anormalement dispersees (defaut de
%                 stabilite process/mesure) -> gonfle MSwithin -> fait chuter
%                 F et ICC sous le seuil ; le retirer les fait repasser au-dessus.
%     'mean'      un specimen a moyenne aberrante (outlier de fabrication) ->
%                 cree/detruit la discrimination selon sa presence.
%
%   opts (struct, optionnel) :
%     .mech    'variance' (defaut) | 'mean'
%     .s,.k    paire cible (defaut s=5, k=3 : descripteur D3 en bord de bande)
%     .i       specimen anormal (defaut 1 : premier specimen de la condition 1)
%     .amp     amplitude de l'anomalie (defaut calibree ; voir run_fragile)
%
%   meta : renvoie s,k,i,mech,amp effectivement utilises (pour tracabilite).

if nargin<3, opts=struct; end
def={'mech','variance';'s',5;'k',3;'i',1;'amp',6.0};
for d=1:size(def,1), if ~isfield(opts,def{d,1}), opts.(def{d,1})=def{d,2}; end, end

% base = dataset nominal (structure et Y identiques au volet 1)
[X, Y] = gen_volet1(P, 'nominal', rs);

G=P.G; n=P.n; R=P.R;
s=opts.s; k=opts.k; i=opts.i; amp=opts.amp;

% moyenne du specimen i a la paire cible (pour inflation centree)
xi = X(i,:,s,k);
mu_i = mean(xi);

switch lower(opts.mech)
    case 'variance'
        % gonfle la dispersion ROI du specimen i, en gardant sa moyenne :
        % chaque ROI est ecarte de sa moyenne d'un facteur amp -> MSwithin explose
        X(i,:,s,k) = mu_i + amp*(xi - mu_i);
    case 'mean'
        % decale la moyenne du specimen i de amp ecarts-types intra
        X(i,:,s,k) = xi + amp*P.sigma_within;
    otherwise
        error('mech inconnu : %s (attendu ''variance'' ou ''mean'')', opts.mech);
end

meta=struct('s',s,'k',k,'i',i,'mech',opts.mech,'amp',amp);
end
