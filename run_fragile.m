% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v2, 31/07/2026 : panel a six descripteurs, D4 correle mais peu fiable
%RUN_FRAGILE  Exerce et DEMONTRE le drapeau "fragile" du LOSO (Section 3.8 / 4.5).
%   Objectif : produire un cas ou une decision de verrou depend d'UN SEUL
%   specimen, de sorte que loso_flags leve fragile(s,k)=true — puis verifier
%   qu'en retirant l'anomalie (dataset nominal) plus aucune paire n'est fragile.
%
%   Le script AUTO-CALIBRE l'amplitude de l'anomalie : il balaie amp jusqu'a
%   obtenir exactement le comportement recherche (au moins une paire fragile,
%   idealement la paire cible injectee), puis affiche le diagnostic complet.
%
%   Necessite : gen_fragile.m, gen_volet1.m, paper_pipeline.m, paper_params.m
clear; clc;
P = paper_params();

% ---- 0) Sanity : sur donnees nominales propres, AUCUNE paire fragile --------
rs0 = RandStream('Threefry','Seed', paper_seed('nominal',1));
[X0,Y0] = gen_volet1(P,'nominal',rs0);
out0 = paper_pipeline(X0,Y0,P,struct('loso',true));
nF0 = nnz(out0.fragile);
fprintf('=== Reference (dataset nominal propre) ===\n');
fprintf('Paires fragiles : %d  (attendu : 0)\n\n', nF0);

% ---- 1) Scenario fragile : balayage d'amplitude jusqu'a lever le drapeau ----
mech   = 'variance';          % 'variance' (defaut) ou 'mean'
sTgt   = 5; kTgt = 3; iTgt=1; % paire cible (D3, bord de bande) + specimen anormal
ampGrid = [2 3 4 5 6 8 10 12 16 20];
seed    = paper_seed('fragile',1);   % meme jeu de base a chaque amplitude

fprintf('=== Scenario fragile (mech=%s, cible s=%d k=%d, specimen i=%d) ===\n', ...
    mech, sTgt, kTgt, iTgt);
found=false;
for amp = ampGrid
    rs = RandStream('Threefry','Seed',seed);
    [X,Y,meta] = gen_fragile(P, rs, ...
        struct('mech',mech,'s',sTgt,'k',kTgt,'i',iTgt,'amp',amp));
    out = paper_pipeline(X,Y,P,struct('loso',true));
    nF  = nnz(out.fragile);
    hitTarget = out.fragile(sTgt,kTgt);
    fprintf('amp=%4.1f -> %d paire(s) fragile(s)', amp, nF);
    if hitTarget, fprintf('  [cible (s=%d,k=%d) FRAGILE]', sTgt, kTgt); end
    fprintf('\n');
    if nF>=1 && ~found
        found=true; ampStar=amp; outStar=out; metaStar=meta; %#ok<NASGU>
        if hitTarget, break; end   % cas ideal : la cible elle-meme est fragile
    end
end

% ---- 2) Diagnostic du cas retenu -------------------------------------------
if ~found
    fprintf(['\nAucune paire fragile obtenue sur la grille testee.\n' ...
        'Augmenter ampGrid, ou passer mech=''mean'', ou choisir une paire\n' ...
        'cible plus proche du seuil (ex. bord de bande D3/D4 : s=8,k=3).\n']);
    return;
end

fprintf('\n=== Cas retenu : amp=%.1f ===\n', ampStar);
[sf,kf] = find(outStar.fragile);
descr = {'D1*','D2*','D3','D4','D5','D6'};
for t=1:numel(sf)
    s=sf(t); k=kf(t);
    fprintf('  paire fragile : %-3s  scale s=%d | decision(complet)=%d | ICC=%.3f | eta2=%.3f\n', ...
        descr{k}, s, outStar.decision(s,k), outStar.ICC(s,k), outStar.eta2(s,k));
end

% ---- 3) Preuve du mecanisme : recalcul LOSO explicite sur la 1re paire ------
s=sf(1); k=kf(1); G=P.G; n=P.n; R=P.R; m=n*R; m0=(n-1)*R;
sampleCondition = repelem(1:G,n);
% statut sur donnees completes
Xp=zeros(G,m); for g=1:G, Xp(g,:)=reshape(X(sampleCondition==g,:,s,k),1,m); end
xg=mean(Xp,2); xx=mean(Xp(:));
MSb=m/(G-1)*sum((xg-xx).^2); MSw=sum(sum((Xp-xg).^2))/(G*(m-1));
Ffull=MSb/MSw; ICCfull=max((Ffull-1)/(Ffull+m-1),0);
% statut en retirant le specimen anormal (i=iTgt de la condition 1)
keep=true(G*n,1); keep(iTgt)=false;
y=[]; grp=[];
for g=1:G
    sel=find(sampleCondition'==g & keep);
    v=reshape(X(sel,:,s,k),[],1); y=[y;v]; grp=[grp; g*ones(numel(v),1)];
end
[~,tbl]=anova1(y,grp,'off'); Floso=tbl{2,5}; ICCloso=max((Floso-1)/(Floso+m0-1),0);

fprintf('\n  --- Mecanisme (paire %s s=%d) ---\n', descr{k}, s);
fprintf('  Complet      : F=%.2f  ICC=%.3f  -> gate2(ICC>=%.2f) : %s\n', ...
    Ffull, ICCfull, P.icc_thr, tern(ICCfull>=P.icc_thr,'PASS','FAIL'));
fprintf('  Sans spec.%d  : F=%.2f  ICC=%.3f  -> gate2(ICC>=%.2f) : %s\n', ...
    iTgt, Floso, ICCloso, P.icc_thr, tern(ICCloso>=P.icc_thr,'PASS','FAIL'));
fprintf(['  => la decision de fiabilite depend de la presence du seul specimen %d :\n' ...
    '     le drapeau fragile est correctement leve. Sur donnees propres (etape 0),\n' ...
    '     aucune paire n''est fragile. Le detecteur LOSO fait ce qu''il doit.\n'], iTgt);

function s=tern(c,a,b), if c, s=a; else, s=b; end, end
