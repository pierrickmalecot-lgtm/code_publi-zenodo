% version v4, 17/09/2026 : couverture des intervalles corriges, analytique et bootstrap
% version v3, 11/09/2026 : champ renomme kappa, seuil KAPPA_MIN a 0.50
% version v2, 11/09/2026 : NREP ramene a 999, borne de paper_seed
% version v1, 05/09/2026 : validation de la correction d'attenuation
%RUN_ATTENUATION  La correction d'attenuation retablit-elle la pente vraie en
%   configuration (b), et jusqu'ou.
%
%   CONTEXTE. En configuration (b), une seule mesure fonctionnelle par specimen,
%   la correlation se fait par moindres carres sur le descripteur moyenne sur
%   les R regions. Cette moyenne etant entachee d'erreur, la pente estimee est
%   tiree vers zero. Le papier signale cette attenuation sans la quantifier.
%   attenuation_correct l'estime et la corrige ; ce script mesure ce que vaut
%   la correction et ou elle cesse d'etre utilisable.
%
%   PROTOCOLE. La reponse est construite sur la valeur VRAIE du specimen, non
%   sur sa moyenne mesuree, de sorte que la pente vraie soit connue et vaille
%   P.beta1. La dispersion entre regions est balayee : plus elle est forte, plus
%   la moyenne est imprecise et plus l'attenuation est severe.
%
%   LECTURE. Deux resumes sont donnes pour la pente corrigee, la moyenne et la
%   mediane. L'ecart entre les deux est le point important : la correction etant
%   un rapport, sa distribution developpe une queue droite lourde des que la
%   fiabilite estimee approche de zero. La mediane reste alors juste tandis que
%   la moyenne diverge. C'est ce qui delimite le domaine d'emploi.
%
%   COUVERTURE. Part des replicats corriges ou l'intervalle contient la pente
%   vraie, pour ci_corr (kappa tenu pour connu) et ci_boot (variance de kappa
%   estimee par bootstrap, BBOOT tirages). Le bootstrap puise dans le flux du
%   replicat APRES la generation des donnees : les jeux simules sont
%   identiques a la version v3.
%
%   Necessite : paper_params, attenuation_correct.
clear; clc;
pool = gcp('nocreate'); if isempty(pool), parpool('Processes'); end

% 999 est le maximum admis par paper_seed : au-dela, le numero de replicat
% empieterait sur l'indice de balayage dans la composition de la graine.
NREP  = 999;
SWGRID = [0.10 0.20 0.30 0.40 0.50 0.60];   % dispersion entre regions
BBOOT  = 1000;                              % tirages bootstrap pour kappa

P = paper_params();
G = P.G; n = P.n; R = P.R; N = P.N;
cond = repelem(P.condLevels(:).', n).';      % niveau vrai de chaque specimen
grp  = repelem((1:G).', n);                  % condition de chaque specimen

nS = numel(SWGRID);
lamM = zeros(nS,1); lamMin = lamM;
bObs = lamM; bCorMean = lamM; bCorMed = lamM; bCorQ95 = lamM;
r2Obs = lamM; r2Cor = lamM; fracOK = lamM;
covA = lamM; covB = lamM; wA = lamM; wB = lamM;

fprintf('=== Correction d attenuation, %d replicats par point ===\n', NREP);
fprintf('Pente vraie beta1 = %.2f, plan G=%d n=%d R=%d\n\n', P.beta1, G, n, R);
tt = tic;

for a = 1:nS
    sw = SWGRID(a);
    bo = zeros(NREP,1); bc = bo; la = bo; ro = bo; rc = bo; ok = bo;
    ca = bo; cb = bo; la2 = bo; lb2 = bo;
    parfor rep = 1:NREP
        rs  = RandStream('Threefry','Seed', paper_seed('nominal', rep, 260+a));
        ai  = P.sigma_between * randn(rs, N, 1);
        tau = cond + ai;                                  % valeur vraie
        X   = tau + sw*randn(rs, N, R);                   % R regions mesurees
        Y   = P.beta0 + P.beta1*tau + P.sigma_func*randn(rs, N, 1);
        o = attenuation_correct(X, Y, 0.05, grp, BBOOT, rs);
        bo(rep) = o.b1; la(rep) = o.kappa; ro(rep) = o.R2;
        ok(rep) = double(o.reliable);
        if o.reliable
            bc(rep) = o.b1_corr; rc(rep) = o.R2_corr;
            ca(rep)  = double(o.ci_corr(1) <= P.beta1 && P.beta1 <= o.ci_corr(2));
            cb(rep)  = double(o.ci_boot(1) <= P.beta1 && P.beta1 <= o.ci_boot(2));
            la2(rep) = diff(o.ci_corr); lb2(rep) = diff(o.ci_boot);
        else
            bc(rep) = NaN; rc(rep) = NaN;
            ca(rep) = NaN; cb(rep) = NaN; la2(rep) = NaN; lb2(rep) = NaN;
        end
    end
    v = ~isnan(bc);
    lamM(a)=median(la); lamMin(a)=min(la); fracOK(a)=mean(ok);
    bObs(a)=median(bo); bCorMean(a)=mean(bc(v)); bCorMed(a)=median(bc(v));
    bCorQ95(a)=quantile(bc(v),0.95);
    r2Obs(a)=median(ro); r2Cor(a)=median(rc(v));
    covA(a)=mean(ca(v)); covB(a)=mean(cb(v));
    wA(a)=median(la2(v)); wB(a)=median(lb2(v));
end
fprintf('Termine en %.1f min.\n\n', toc(tt)/60);

% --- sorties ------------------------------------------------------------------
fprintf('--- Pente : attenuation observee et correction ---\n');
fprintf('%8s %9s %10s %10s %12s %12s %10s\n', ...
    'sigma_w','kappa','attenuat.','pente obs','corr. med','corr. moy','corr. q95');
for a = 1:nS
    fprintf('%8.2f %9.3f %9.1f %% %10.3f %12.3f %12.3f %10.2f\n', ...
        SWGRID(a), lamM(a), 100*(1-lamM(a)), bObs(a), ...
        bCorMed(a), bCorMean(a), bCorQ95(a));
end

fprintf('\n--- Biais relatif a la pente vraie ---\n');
fprintf('%8s %14s %14s %14s\n','sigma_w','non corrigee','corrigee med','corrigee moy');
for a = 1:nS
    fprintf('%8.2f %13.1f %% %13.1f %% %13.1f %%\n', SWGRID(a), ...
        100*(bObs(a)/P.beta1-1), 100*(bCorMed(a)/P.beta1-1), ...
        100*(bCorMean(a)/P.beta1-1));
end

fprintf('\n--- Coefficient de determination et applicabilite ---\n');
fprintf('%8s %12s %12s %16s %12s\n', ...
    'sigma_w','R2 observe','R2 corrige','kappa minimal','corrigible');
for a = 1:nS
    fprintf('%8.2f %12.3f %12.3f %16.3f %11.1f %%\n', SWGRID(a), ...
        r2Obs(a), r2Cor(a), lamMin(a), 100*fracOK(a));
end

fprintf('\n--- Couverture a 95 %% des intervalles corriges ---\n');
fprintf('%8s %14s %14s %14s %14s\n', ...
    'sigma_w','couv. analyt.','couv. boot','larg. analyt.','larg. boot');
for a = 1:nS
    fprintf('%8.2f %13.1f %% %13.1f %% %14.2f %14.2f\n', SWGRID(a), ...
        100*covA(a), 100*covB(a), wA(a), wB(a));
end

fprintf(['\nLECTURE. La colonne corrigee med doit rester proche de zero de biais\n' ...
         'tant que l attenuation est moderee. L ecart entre corrigee med et\n' ...
         'corrigee moy mesure la lourdeur de la queue droite : au-dela d une\n' ...
         'certaine dispersion entre regions, la moyenne diverge alors que la\n' ...
         'mediane tient, et la correction ne doit plus etre rapportee sans\n' ...
         'reserve. La colonne corrigible donne la part de replicats ou la\n' ...
         'fiabilite estimee reste au-dessus du seuil.\n']);

save('ws_run_attenuation.mat','SWGRID','lamM','lamMin','bObs','bCorMean', ...
     'bCorMed','bCorQ95','r2Obs','r2Cor','fracOK','covA','covB','wA','wB', ...
     'BBOOT','NREP','P');
fprintf('\nws_run_attenuation.mat ecrit.\n');
