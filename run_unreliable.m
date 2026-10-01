% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 31/07/2026 : scenario dedie du descripteur correle mais peu fiable
%RUN_UNRELIABLE  Scenario dedie : un descripteur correle mais peu fiable.
%
%   POURQUOI UN SCENARIO SEPARE. Le panel de reference ne contient aucun
%   descripteur qui soit associe a la reponse tout en etant mesure trop
%   imprecisement pour fonder une specification, c'est-a-dire exactement la
%   situation que le verrou de fiabilite existe pour intercepter. Produire ce
%   cas demande un terme local partage entre le descripteur et la reponse, sans
%   contrepartie physique mesurable sur une surface reelle. Le laisser dans le
%   panel principal aurait abaisse tous les R2 du jeu de reference pour une
%   construction difficile a defendre. Il est donc traite ici a part, au meme
%   titre que le tirage independant du scenario 'null'.
%
%   CONSTRUCTION. Un sixieme descripteur est ajoute au panel. Sa moyenne suit la
%   condition procede sur les echelles 5 a 7, comme D1, et il porte en plus un
%   terme z tire une fois par region de mesure, d'ecart-type sigma_z. Le meme z
%   entre dans la reponse avec le coefficient gamma_z. Le descripteur et la
%   reponse partagent donc une cause commune sans que l'un agisse sur l'autre.
%   Ce terme gonfle la dispersion entre regions de mesure, donc fait chuter la
%   fiabilite, tout en etant predictif de la reponse, donc en soutenant la
%   correlation. Aucun autre mecanisme n'y parvient : toute variabilite propre
%   au descripteur fait chuter le R2 plus vite que l'ICC, la pente s'attenuant
%   et le residu grossissant en meme temps.
%
%   MESURES. Le R2 du sixieme descripteur est calcule deux fois, une fois par la
%   procedure complete, ou il n'existe pas puisque la paire est arretee au
%   verrou de fiabilite, et une fois avec les verrous amont neutralises, ce qui
%   donne la valeur que la paire aurait obtenue si rien ne l'avait arretee.
%   C'est cette seconde valeur qui etablit qu'elle aurait ete retenue a tort.

clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

P = paper_params();
P.K = 6;                      % le scenario ajoute un descripteur au panel
kU  = 6; band = 5:7;          % le descripteur correle mais peu fiable
names = {'D1','D2','D3','D4','D5','DU'};
NREP = 200;

% --- tirage de reference ------------------------------------------------------
rs = RandStream('Threefry','Seed',40007);
[X,Y] = gen_volet1(P,'unreliable',rs);
out = paper_pipeline(X,Y,P,struct('bootstrap',true,'B',1000));

% R2 que la paire aurait obtenu sans les verrous amont
Pn = P; Pn.q=1; Pn.eta2_thr=0; Pn.icc_thr=0; Pn.r2_thr=-Inf;
outF = paper_pipeline(X,Y,Pn);

fprintf('=== Descripteur correle mais peu fiable, tirage de reference ===\n\n');
fprintf('%-6s %8s %8s %8s %10s %9s %10s\n', ...
    'scale','F','eta2','ICC','R2 force','rho','decision');
for s = 1:P.S
    fprintf('s=%-4d %8.1f %8.3f %8.3f %10.3f %9.3f %10d\n', s, ...
        out.F(s,kU), out.eta2(s,kU), out.ICC(s,kU), outF.R2m(s,kU), ...
        out.rho(s,kU), out.decision(s,kU));
end

fprintf('\n=== Panel complet, tirage de reference ===\n');
for k = 1:P.K
    fprintf('  %-4s retenu sur %s\n', names{k}, mat2str(find(out.decision(:,k)==3)'));
end

% --- stabilite sur replicats --------------------------------------------------
dec = zeros(P.S,4); r2f = zeros(NREP,1); ret = zeros(NREP,1);
iccU = zeros(NREP,1); rhoU = zeros(NREP,1);
parfor r = 1:NREP
    rr = RandStream('Threefry','Seed', paper_seed('unreliable', r));
    [Xr,Yr] = gen_volet1(P,'unreliable',rr);
    o  = paper_pipeline(Xr,Yr,P);
    Pnn = P; Pnn.q=1; Pnn.eta2_thr=0; Pnn.icc_thr=0; Pnn.r2_thr=-Inf;
    of = paper_pipeline(Xr,Yr,Pnn);
    loc = zeros(P.S,4);
    for s = 1:P.S, loc(s, o.decision(s,kU)+1) = 1; end
    dec  = dec + loc;
    r2f(r)  = mean(of.R2m(band,kU),'omitnan');
    iccU(r) = mean(o.ICC(band,kU));
    rhoU(r) = mean(o.rho(band,kU));
    ret(r)  = nnz(o.decision==3);
end
dec = dec / NREP;

fprintf('\n=== Stabilite sur %d replicats ===\n', NREP);
fprintf('%-6s %10s %10s %10s %10s\n','scale','exclu','discrim.','+fiable','retenu');
for s = 1:P.S
    fprintf('s=%-4d %10.2f %10.2f %10.2f %10.2f\n', s, dec(s,:));
end

fprintf('\nSur les echelles %d a %d :\n', band(1), band(end));
fprintf('  ICC moyen                     %.3f\n', mean(iccU));
fprintf('  rho moyen                     %.3f\n', mean(rhoU));
fprintf('  R2 force, moyenne             %.3f\n', mean(r2f));
fprintf('  R2 force au-dessus de %.2f     %.0f pour cent des replicats\n', ...
    P.r2_thr, 100*mean(r2f > P.r2_thr));
fprintf('  paires retenues par replicat, panel entier : %.2f\n', mean(ret));

fprintf(['\nLecture. Le descripteur est discriminant sur ses echelles et le R2\n' ...
    'qu''il aurait obtenu depasse le seuil de correlation, donc la procedure\n' ...
    'sans verrou amont l''aurait retenu. Son ICC le place pourtant loin sous\n' ...
    'le seuil de fiabilite, et la procedure complete l''ecarte. C''est la\n' ...
    'seule configuration ou le verrou de fiabilite est indispensable : il\n' ...
    'intercepte une association reelle mais non specifiable.\n']);

save('unreliable_results.mat','dec','r2f','iccU','rhoU','ret','out','outF','NREP');
fprintf('\nResultats enregistres dans unreliable_results.mat\n');
