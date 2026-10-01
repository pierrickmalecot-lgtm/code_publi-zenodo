% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v2, 04/09/2026 : ajout du scenario null, le vrai test du mode d echec
% version v1, 04/09/2026 : validation du drapeau de non-monotonie
%RUN_USHAPE  Le drapeau de non-monotonie se leve-t-il quand il doit, et
%   reste-t-il muet sinon.
%
%   LE PROBLEME. Une relation a optimum interieur, en U ou en cloche, est
%   courante en tribologie : le frottement peut etre minimal a rugosite
%   intermediaire, l'adherence maximale a une valeur intermediaire de portance.
%   Sur une telle paire, la pente lineaire s'annule et le Spearman aussi, donc
%   la porte 3 rejette. Rien dans la sortie n'indique que la relation existe
%   mais n'est pas de la forme ajustee : la paire est ecartee en silence.
%
%   LE SIGNALEMENT. Toutes les quantites necessaires sont deja calculees. Une
%   paire ayant franchi les portes 1 et 2 discrimine les conditions et est
%   mesuree de facon reproductible : une structure existe. Si son R2 marginal
%   est neanmoins sous le seuil ET que son Spearman est faible en valeur
%   absolue, c'est la FORME ajustee qui est en cause, pas l'absence de
%   relation. Le drapeau est exactement cette conjonction. Il ne modifie
%   aucune decision, il ne s'applique qu'aux paires deja rejetees.
%
%   CE QUE CE SCRIPT ETABLIT.
%     (1) SENSIBILITE. Sur le scenario 'ushape' a optimum centre, le drapeau
%         se leve sur la paire porteuse (D1, s = s_star).
%     (2) SPECIFICITE, deux cas a distinguer.
%         'nominal' et 'curvature' : la paire FRANCHIT la porte 3, donc la
%         condition ~gate3 l'exclut d'emblee et le drapeau ne peut pas se
%         lever. Ces deux lignes ne testent pas le seuil sur le Spearman.
%         'null' : la reponse est tiree independamment de la topographie. Les
%         descripteurs discriminants franchissent les portes 1 et 2, echouent
%         a la porte 3, et leur Spearman est proche de zero. Ils remplissent
%         donc la condition du drapeau SANS qu'aucune relation n'existe.
%         C'est le vrai test du mode d'echec, et aucune valeur du seuil ne
%         separe ce cas d'une cloche centree : les deux donnent R2 bas et
%         Spearman bas.
%     (3) DOMAINE. En deplacant l'optimum vers un bord de la plage, une
%         composante monotone reapparait, le Spearman remonte et le drapeau
%         retombe. Le balayage situe la transition.
%     (4) REMEDE ET TRI. Sur les paires signalees, le reajustement quadratique
%         de la section 3.9.6 est ce qui separe une vraie forme non monotone
%         d'une absence de relation : le drapeau ratisse large, l'ecart d'AIC
%         tranche. C'est pourquoi il est mesure sur 'ushape' ET sur 'null'.
%
%   Necessite : gen_volet1 v8, paper_pipeline v9, paper_params v8.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP  = 200;
X0    = [0.50 0.48 0.46 0.44 0.42 0.40 0.35 0.30];  % position de l'optimum

P0 = paper_params();
S = P0.S; K = P0.K; ks = 1; ss = P0.s_star;   % paire porteuse : D1 au pic

% --- (1) et (2) sensibilite et specificite ------------------------------------
SCEN = {'ushape','nominal','curvature','null'};
CURV = [0 0 3 0];               % taux de courbure, actif seulement sur le 3e
flag_star = zeros(numel(SCEN),1);   % drapeau sur la paire porteuse
flag_tot  = zeros(numel(SCEN),1);   % drapeaux sur toute la grille
r2_star   = zeros(numel(SCEN),1);
rho_star  = zeros(numel(SCEN),1);
ret_tot   = zeros(numel(SCEN),1);
dAIC_med  = nan(numel(SCEN),1);     % ecart d'AIC median sur les paires signalees
dAIC_pos  = nan(numel(SCEN),1);     % fraction des signalees ou l'AIC favorise
                                    % le quadratique

fprintf('=== Drapeau de non-monotonie, %d replicats ===\n\n', NREP);
tt = tic;
for is = 1:numel(SCEN)
    P = P0; P.curv = CURV(is); P.ushape_x0 = 0.50;
    fs = zeros(NREP,1); ft = fs; r2 = fs; rr = fs; rt = fs;
    da = cell(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(SCEN{is}, r));
        [X,Y] = gen_volet1(P, SCEN{is}, rs);
        o = paper_pipeline(X,Y,P,struct('quadCheck',true));
        fs(r) = double(o.nonmono(ss,ks));
        ft(r) = sum(o.nonmono(:));
        r2(r) = o.R2m(ss,ks);
        rr(r) = abs(o.rho(ss,ks));
        rt(r) = sum(o.decision(:)==3);
        v = o.dAIC(o.nonmono);          % toutes les paires signalees
        da{r} = v(~isnan(v));
    end
    flag_star(is)=mean(fs); flag_tot(is)=mean(ft);
    r2_star(is)=mean(r2,'omitnan'); rho_star(is)=mean(rr);
    ret_tot(is)=mean(rt);
    allda = vertcat(da{:});
    if ~isempty(allda)
        dAIC_med(is) = median(allda);
        dAIC_pos(is) = mean(allda > 0);
    end
    fprintf('%-12s termine (%.1f min)\n', SCEN{is}, toc(tt)/60);
end

fprintf('\n--- Sensibilite et specificite, paire porteuse D1 a s = %d ---\n', ss);
fprintf('%-12s %9s %8s %11s %11s %9s %10s %9s\n', ...
    'scenario','R2 marg','|rho|','drapeau D1*','drapeaux/50','retenues', ...
    'dAIC med','AIC quad');
for is = 1:numel(SCEN)
    fprintf('%-12s %9.3f %8.3f %11.3f %11.2f %9.2f %10.1f %8.1f %%\n', SCEN{is}, ...
        r2_star(is), rho_star(is), flag_star(is), flag_tot(is), ret_tot(is), ...
        dAIC_med(is), 100*dAIC_pos(is));
end
fprintf(['\nLECTURE. Le drapeau doit valoir ~1 sur ushape. Sur nominal et\n' ...
         'curvature il ne peut pas se lever, la paire franchissant la porte 3.\n' ...
         'Sur null il DOIT se lever : le drapeau seul ne distingue pas une\n' ...
         'cloche centree d une absence de relation. C est l ecart d AIC qui\n' ...
         'les separe, positif sur ushape et negatif ou nul sur null.\n' ...
         'La colonne retenues doit rester celle de chaque scenario : le\n' ...
         'drapeau ne retient ni n exclut aucune paire.\n']);

% --- (3) domaine : balayage de la position de l'optimum ------------------------
fprintf('\n=== Balayage de la position de l optimum ===\n');
fl = zeros(numel(X0),1); r2b = fl; rob = fl; retb = fl;
for ix = 1:numel(X0)
    P = P0; P.ushape_x0 = X0(ix);
    f = zeros(NREP,1); a = f; b = f; c = f;
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('ushape', r, 400+ix));   % bloc 400-419
        [X,Y] = gen_volet1(P,'ushape',rs);
        o = paper_pipeline(X,Y,P);
        f(r) = double(o.nonmono(ss,ks));
        a(r) = o.R2m(ss,ks);
        b(r) = abs(o.rho(ss,ks));
        c(r) = double(o.decision(ss,ks)==3);
    end
    fl(ix)=mean(f); r2b(ix)=mean(a,'omitnan'); rob(ix)=mean(b); retb(ix)=mean(c);
end

fprintf('%8s %10s %10s %12s %14s\n', ...
    'x0','R2 marg','|rho|','drapeau','paire retenue');
for ix = 1:numel(X0)
    fprintf('%8.2f %10.3f %10.3f %12.3f %14.3f\n', ...
        X0(ix), r2b(ix), rob(ix), fl(ix), retb(ix));
end
fprintf(['\nLECTURE. Aux optima centres, les deux indices sont bas et le drapeau\n' ...
         'se leve. En deplacant l optimum vers un bord, une composante monotone\n' ...
         'reapparait, le Spearman remonte, le drapeau retombe et la paire finit\n' ...
         'par etre retenue par la porte 3 elle-meme.\n']);

% --- (4) remede et tri : ushape contre null -----------------------------------
% Le drapeau seul ne separe pas une cloche centree d'une absence de relation.
% Ce volet mesure ce que le reajustement quadratique ajoute sur les deux, en
% suivant la paire porteuse D1 a s = s_star dans les deux scenarios.
fprintf('\n=== Reajustement quadratique : ushape contre null ===\n');
QS = {'ushape','null'};
fprintf('%-10s %10s %12s %12s %12s %11s\n', ...
    'scenario','signalee','R2 lineaire','R2 quadrat.','dAIC median','AIC quad');
qres = struct();
for iq = 1:numel(QS)
    P = P0; P.ushape_x0 = 0.50;
    q_r2 = nan(NREP,1); q_da = nan(NREP,1); q_lin = nan(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(QS{iq}, r, 420+iq));   % bloc 420-439
        [X,Y] = gen_volet1(P, QS{iq}, rs);
        o = paper_pipeline(X,Y,P,struct('quadCheck',true));
        if o.nonmono(ss,ks)
            q_lin(r) = o.R2m(ss,ks);
            q_r2(r)  = o.R2m_quad(ss,ks);
            q_da(r)  = o.dAIC(ss,ks);
        end
    end
    v = ~isnan(q_da);
    fprintf('%-10s %9.1f %% %12.3f %12.3f %12.1f %9.1f %%\n', QS{iq}, ...
        100*mean(v), mean(q_lin(v)), mean(q_r2(v)), median(q_da(v)), ...
        100*mean(q_da(v) > 0));
    qres.(QS{iq}) = struct('flagged',mean(v),'r2lin',mean(q_lin(v)), ...
        'r2quad',mean(q_r2(v)),'dAIC',q_da(v));
end
fprintf(['\nLECTURE. Le drapeau se leve dans les deux cas. Le reajustement les\n' ...
         'separe : sur ushape le quadratique recupere la relation et l AIC le\n' ...
         'designe, sur null il n explique rien de plus et l AIC le rejette.\n' ...
         'C est donc l ecart d AIC, non le drapeau, qui porte le diagnostic.\n']);

save('ws_run_ushape.mat','flag_star','flag_tot','r2_star','rho_star','ret_tot', ...
     'dAIC_med','dAIC_pos','SCEN','CURV','X0','fl','r2b','rob','retb', ...
     'qres','NREP','P0');
fprintf('\nws_run_ushape.mat ecrit. Total %.1f min.\n', toc(tt)/60);
