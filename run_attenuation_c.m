% version v2, 17/09/2026 : intervalle avec variance bootstrap de kappa
% version v1, 17/09/2026 : validation de la correction d'attenuation, configuration (c)
%RUN_ATTENUATION_C  La correction d'attenuation validee en configuration (b)
%   vaut-elle aussi en configuration (c).
%
%   CONTEXTE. En configuration (c), T essais fonctionnels par specimen, sans
%   appariement aux regions. Le modele mixte garde l'intercept aleatoire
%   specimen et entre le descripteur au niveau specimen, par sa moyenne sur
%   les R regions :
%       Y_it = beta0 + beta1 * Xbar_i + u_i + eps_it
%   Le predicteur est le meme qu'en (b), donc entache de la meme erreur de
%   moyennage. La question est de savoir si la division par kappa, validee en
%   (b) par run_attenuation, retablit aussi la pente du modele mixte.
%
%   ATTENDU. Quand T est le meme pour tous les specimens, la pente du modele
%   mixte a intercept aleatoire, pour un predicteur de niveau specimen, est
%   exactement celle des moindres carres de Ybar_i sur Xbar_i. La configuration
%   (c) se ramene alors a la (b), avec une reponse moyennee sur T essais.
%   L'attenuation ne depend que de l'erreur sur Xbar : kappa est le meme, et la
%   correction devrait s'appliquer telle quelle. Quand T varie d'un specimen a
%   l'autre, l'equivalence n'est plus exacte ; c'est ce cas que le script
%   mesure vraiment. La colonne ecart LME-MCO controle l'equivalence.
%
%   PROTOCOLE. Generateur identique a run_attenuation pour le descripteur. La
%   reponse est construite sur la valeur VRAIE du specimen, avec un effet
%   specimen u_i (P.sigma_u) et un bruit d'essai eps_it (P.sigma_func). Deux
%   plans :
%     equilibre      T = 3 pour tous les specimens
%     desequilibre   T tire uniformement dans {2, 3, 4} pour chaque specimen
%   kappa est estime sur X seul, par attenuation_correct, avec le meme seuil.
%   attenuation_correct n'est pas modifiee : elle recoit Ybar pour fournir
%   kappa, et l'ajustement mixte est fait ici.
%
%   COUVERTURE. Part des replicats corriges ou l'intervalle contient la pente
%   vraie. Deux intervalles, construits sur la pente et l'ecart-type du modele
%   mixte :
%     analytique  b1_corr +/- t * se/kappa, kappa tenu pour connu
%     bootstrap   b1_corr +/- t * sqrt(se^2/kappa^2 + b1^2*var_boot(kappa)/kappa^4)
%   La variance bootstrap de kappa vient de attenuation_correct (BBOOT tirages,
%   stratifies par condition). Le bootstrap puise dans le flux du replicat
%   APRES la generation des donnees : les jeux simules sont identiques a la v1.
%
%   GRAINES. Bloc 440-459 de paper_seed : 440+a pour le plan equilibre,
%   450+a pour le plan desequilibre. Bloc libre au 17/09/2026, a declarer
%   dans l'en-tete de paper_seed.
%
%   Necessite : paper_params, paper_seed, attenuation_correct,
%   Statistics and Machine Learning Toolbox (fitlme).
clear; clc;
pool = gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP   = 999;                               % maximum admis par paper_seed
SWGRID = [0.10 0.20 0.30 0.40 0.50 0.60];   % dispersion entre regions
T_EQ   = 3;                                 % essais par specimen, plan equilibre
T_RNG  = [2 4];                             % bornes de T, plan desequilibre
PLANS  = {'equilibre', 'desequilibre'};
LVL0   = [440 450];
ALPHA  = 0.05;
BBOOT  = 1000;                              % tirages bootstrap pour kappa

P = paper_params();
G = P.G; n = P.n; R = P.R; N = P.N;
cond = repelem(P.condLevels(:).', n).';     % niveau vrai de chaque specimen
grp  = repelem((1:G).', n);                 % condition de chaque specimen
tq   = tinv(1-ALPHA/2, N-2);

nS = numel(SWGRID); nP = numel(PLANS);
kapM = zeros(nS,nP); bObs = kapM; bCorMed = kapM; bCorMean = kapM;
bCorQ95 = kapM; fracOK = kapM; cover = kapM; dEquiv = kapM;
coverB = kapM; wA = kapM; wB = kapM;

tt = tic;
for p = 1:nP
    for a = 1:nS
        sw = SWGRID(a);
        bo = zeros(NREP,1); bc = bo; ka = bo; ok = bo; cv = bo; de = bo;
        cvb = bo; wa = bo; wb = bo;
        parfor rep = 1:NREP
            rs  = RandStream('Threefry','Seed', ...
                             paper_seed('nominal', rep, LVL0(p)+a));
            ai  = P.sigma_between * randn(rs, N, 1);
            tau = cond + ai;                              % valeur vraie
            X   = tau + sw*randn(rs, N, R);               % R regions mesurees
            if p == 1
                Ti = T_EQ * ones(N,1);
            else
                Ti = randi(rs, T_RNG, N, 1);
            end
            ui  = P.sigma_u * randn(rs, N, 1);
            spec = repelem((1:N).', Ti);
            Y   = P.beta0 + P.beta1*tau(spec) + ui(spec) ...
                  + P.sigma_func*randn(rs, sum(Ti), 1);
            xb  = mean(X, 2);

            % ajustement mixte, descripteur au niveau specimen
            tbl = table(Y, xb(spec), categorical(spec), ...
                        'VariableNames', {'Y','Xb','spec'});
            lme = fitlme(tbl, 'Y ~ 1 + Xb + (1|spec)');
            b1  = lme.Coefficients.Estimate(2);
            se  = lme.Coefficients.SE(2);

            % kappa sur X seul ; b1 de attenuation_correct = MCO sur Ybar
            Ybar = accumarray(spec, Y) ./ Ti;
            o = attenuation_correct(X, Ybar, ALPHA, grp, BBOOT, rs);

            bo(rep) = b1; ka(rep) = o.kappa; ok(rep) = double(o.reliable);
            de(rep) = abs(b1 - o.b1);
            if o.reliable
                bcr = b1 / o.kappa;
                bc(rep) = bcr;
                hA  = tq * se / o.kappa;
                hB  = tq * sqrt(se^2/o.kappa^2 + b1^2*o.var_kappa_boot/o.kappa^4);
                cv(rep)  = double(abs(bcr - P.beta1) <= hA);
                cvb(rep) = double(abs(bcr - P.beta1) <= hB);
                wa(rep) = 2*hA; wb(rep) = 2*hB;
            else
                bc(rep) = NaN; cv(rep) = NaN;
                cvb(rep) = NaN; wa(rep) = NaN; wb(rep) = NaN;
            end
        end
        v = ~isnan(bc);
        kapM(a,p) = median(ka); fracOK(a,p) = mean(ok);
        bObs(a,p) = median(bo); dEquiv(a,p) = median(de);
        bCorMed(a,p) = median(bc(v)); bCorMean(a,p) = mean(bc(v));
        bCorQ95(a,p) = quantile(bc(v), 0.95);
        cover(a,p) = mean(cv(v));
        coverB(a,p) = mean(cvb(v));
        wA(a,p) = median(wa(v)); wB(a,p) = median(wb(v));
    end
end
fprintf('Termine en %.1f min.\n\n', toc(tt)/60);

% --- sorties ------------------------------------------------------------------
fprintf('=== Correction d attenuation, configuration (c), %d replicats par point ===\n', NREP);
fprintf('Pente vraie beta1 = %.2f, plan G=%d n=%d R=%d\n', P.beta1, G, n, R);
fprintf('sigma_u = %.2f, sigma_eps = %.2f\n', P.sigma_u, P.sigma_func);

for p = 1:nP
    if p == 1
        fprintf('\n##### Plan %s, T = %d #####\n', PLANS{p}, T_EQ);
    else
        fprintf('\n##### Plan %s, T dans {%d..%d} #####\n', PLANS{p}, T_RNG);
    end

    fprintf('\n--- Pente : attenuation observee et correction ---\n');
    fprintf('%8s %9s %10s %10s %12s %12s %10s\n', ...
        'sigma_w','kappa','attenuat.','pente obs','corr. med','corr. moy','corr. q95');
    for a = 1:nS
        fprintf('%8.2f %9.3f %9.1f %% %10.3f %12.3f %12.3f %10.2f\n', ...
            SWGRID(a), kapM(a,p), 100*(1-kapM(a,p)), bObs(a,p), ...
            bCorMed(a,p), bCorMean(a,p), bCorQ95(a,p));
    end

    fprintf('\n--- Biais relatif a la pente vraie ---\n');
    fprintf('%8s %14s %14s %14s\n','sigma_w','non corrigee','corrigee med','corrigee moy');
    for a = 1:nS
        fprintf('%8.2f %13.1f %% %13.1f %% %13.1f %%\n', SWGRID(a), ...
            100*(bObs(a,p)/P.beta1-1), 100*(bCorMed(a,p)/P.beta1-1), ...
            100*(bCorMean(a,p)/P.beta1-1));
    end

    fprintf('\n--- Applicabilite et equivalence ---\n');
    fprintf('%8s %12s %16s\n','sigma_w','corrigible','ecart LME-MCO');
    for a = 1:nS
        fprintf('%8.2f %11.1f %% %16.2e\n', SWGRID(a), ...
            100*fracOK(a,p), dEquiv(a,p));
    end

    fprintf('\n--- Couverture a 95 %% des intervalles corriges ---\n');
    fprintf('%8s %14s %14s %14s %14s\n', ...
        'sigma_w','couv. analyt.','couv. boot','larg. analyt.','larg. boot');
    for a = 1:nS
        fprintf('%8.2f %13.1f %% %13.1f %% %14.2f %14.2f\n', SWGRID(a), ...
            100*cover(a,p), 100*coverB(a,p), wA(a,p), wB(a,p));
    end
end

fprintf(['\nLECTURE. Plan equilibre : l ecart LME-MCO doit etre de l ordre de la\n' ...
         'precision numerique, et les tableaux proches de ceux de run_attenuation.\n' ...
         'Plan desequilibre : c est le cas a juger. Si le biais corrige en mediane\n' ...
         'reste du meme ordre qu en (b), la correction vaut pour (c). La couverture\n' ...
         'bootstrap doit etre proche de 95 %%. L analytique, qui tient kappa pour\n' ...
         'connu, est attendu trop etroit des que l attenuation depasse 20 %%.\n']);

save('ws_run_attenuation_c.mat','SWGRID','PLANS','T_EQ','T_RNG','kapM', ...
     'bObs','bCorMed','bCorMean','bCorQ95','fracOK','cover','coverB','wA','wB', ...
     'dEquiv','BBOOT','NREP','P');
fprintf('\nws_run_attenuation_c.mat ecrit.\n');
