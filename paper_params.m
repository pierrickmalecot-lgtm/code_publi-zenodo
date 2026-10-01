% version v8, 04/09/2026 : scenario 'ushape' et drapeau de non-monotonie
% version v7, 02/08/2026 : ajout de rho_spat, curv et sigma_u1
% version v6, 31/07/2026 : panel de reference a cinq descripteurs
function P = paper_params()
%PAPER_PARAMS  Parametres du volet 1 (papier, section 4.1).
%   VERSION 6 : le panel de reference compte CINQ descripteurs. Le second
%   discriminant seul a ete supprime, redondant avec D3 dont il ne differait
%   que par la bande visee. Le descripteur correle mais peu fiable, qui
%   l'avait un temps remplace, fait desormais l'objet d'un scenario dedie :
%   le terme local qui le construit n'a pas de contrepartie physique
%   mesurable, et le laisser dans le panel de reference abaissait tous les R2
%   du jeu principal sans raison defendable.
P.G=5; P.n=5; P.R=5; P.S=10; P.K=5; P.N=P.G*P.n;
P.condLevels = linspace(0.2, 1.0, P.G);
P.sigma_between = 0.15;   % a_i
P.sigma_within  = 0.10;   % w
P.sigma_func    = 0.20;   % eps (Y)
P.sigma_u       = 0.25;   % intercept aleatoire
P.beta0=2.0; P.beta1=3.5; P.s_star=6;
% Reglages du scenario 'unreliable', ou un sixieme descripteur est ajoute au
% panel. Un terme local z, propre a la region de mesure, lui est ajoute et
% entre aussi dans la reponse avec le coefficient gamma_z. sigma_z pilote la
% fiabilite, gamma_z la correlation. Ces valeurs n'ont aucun effet sur les
% autres scenarios, ou le terme local est nul.
% Reglage verifie sur 200 replicats : ICC 0.22, R2 force 0.53 en moyenne,
% au-dessus du seuil de correlation dans la totalite des replicats.
P.sigma_z       = 0.55;
P.gamma_z       = 1.50;
% Parametres des tests de robustesse (section 6.4). Les deux valeurs par
% defaut sont neutres : rho_spat = 0 rend le bruit AR(1) identique au bruit
% independant, curv = 0 rend la reponse saturante identique a la
% reponse lineaire. Elles n'ont donc aucun effet sur les scenarios existants
% et ne modifient aucun resultat deja publie.
P.rho_spat = 0.0;   % autocorrelation AR(1) du bruit ROI le long de r
P.curv     = 0.0;   % taux de courbure de la reponse saturante
P.sigma_u1 = 0.0;   % ecart-type de la pente aleatoire (scenario 'randslope')
% Scenario 'ushape' : reponse a optimum interieur, cas des relations en U ou
% en cloche courantes en tribologie (frottement minimal a rugosite
% intermediaire). x0 situe la position de l'optimum dans la plage normalisee
% [0,1] du descripteur porteur. A x0 = 0.5 la relation est symetrique : la
% pente lineaire et le Spearman sont nuls tous les deux, et la paire est
% rejetee a la porte 3 sans que rien ne signale pourquoi. En deplacant x0 vers
% un bord, une composante monotone reapparait et le Spearman remonte. Ce
% parametre n'a d'effet que sous le scenario 'ushape'.
P.ushape_x0 = 0.50;
% Scenario 'latent' : ecart-type du bruit separant le descripteur porteur de la
% variable latente qui produit la reponse, en fraction de l'ecart-type de X sur
% la paire porteuse. A 0 le scenario est identique a 'nominal'. Ce parametre
% n'a d'effet que sous ce scenario.
P.lat_noise = 0.50;
% Drapeau de non-monotonie (inspection, pas une porte). Une paire ayant
% franchi les deux premieres portes, dont le R2 marginal est sous le seuil de
% correlation ET dont le Spearman est faible en valeur absolue, est le
% signalement d'une relation non monotone : la porte 3 la rejette alors que la
% relation existe, mais n'est pas de la forme ajustee. Voir section 3.9.6.
P.rho_nonmono = 0.30;
% Seuils (papier)
P.q=0.05; P.eta2_thr=0.06; P.icc_thr=0.60; P.r2_thr=0.30; P.C_thr=0.90;
end
