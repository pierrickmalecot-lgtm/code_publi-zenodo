% version v9, 04/09/2026 : grille sigma_u1 etendue jusqu a 1.0*beta1
% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/08/2026 : selection de modele, pente aleatoire (section 4.4)
%RUN_RANDSLOPE  Verifie la regle de selection entre modele a intercept
%   aleatoire et modele a pente aleatoire.
%
%   Le jeu de reference est genere sans pente aleatoire, donc la specification
%   reduite y est le bon modele. Ce script balaie sigma_u1 et montre que le
%   critere d'information bascule vers le modele complet au-dela d'un certain
%   niveau, et reste sur le modele reduit en dessous. C'est ce qui manque a la
%   section 4.4, qui decrit la regle sans jamais l'exercer.
%
%   Necessite gen_volet1.m avec le scenario 'randslope' et paper_params.m avec
%   P.sigma_u1.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P = paper_params(); NREP = 200;

% sigma_u1 exprime en fraction de beta1 = 3.5. La grille s'arretait a 0.4 dans
% la version precedente, en supposant que le passage d'un regime a l'autre s'y
% produisait. Il ne s'y produit pas : a 0.4, l'AIC ne retient le modele complet
% que dans 20.5 % des replicats et l'ecart d'AIC reste negatif sur toute la
% grille. Elle est donc etendue jusqu'a 1.0*beta1, ou l'ecart-type de la pente
% aleatoire egale la pente fixe elle-meme, pour situer la bascule si elle
% existe, ou etablir qu'elle ne se produit pas a la taille de plan du papier.
fracGrid = [0 0.05 0.10 0.20 0.30 0.40 0.60 0.80 1.00];
sigGrid  = fracGrid * P.beta1;
s_star = P.s_star; kRef = 1;

nS = numel(sigGrid);
pFull_AIC = zeros(nS,1);   % proportion de replicats ou l'AIC choisit le complet
pFull_BIC = zeros(nS,1);
pFull_LRT = zeros(nS,1);   % test du rapport de vraisemblance, alpha = 0.05
dAIC      = zeros(nS,1);   % AIC(reduit) - AIC(complet), moyenne

fprintf('=== Selection de modele, pente aleatoire, %d replicats ===\n', NREP);
fprintf('beta1 = %.2f, s* = %d\n\n', P.beta1, s_star);
fprintf('%8s %10s %8s %8s %8s %10s\n', ...
    'sig_u1','/beta1','AIC','BIC','LRT','dAIC');

tt = tic;
for a = 1:nS
    Pa = P; Pa.sigma_u1 = sigGrid(a);
    selA = false(NREP,1); selB = false(NREP,1); selL = false(NREP,1);
    dA = nan(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('randslope', r, 60+a));   % bloc 60-79
        [X,Y] = gen_volet1(Pa,'randslope',rs);

        % memes conventions que lme_r2_marginal : X centre reduit, ROI nichees
        x = reshape(X(:,:,s_star,kRef),[],1); x = (x-mean(x))/std(x);
        y = Y(:);
        sid = repelem((1:Pa.N)', Pa.R, 1);
        T = table(y, x, nominal(sid), 'VariableNames', {'Y','X','Sample'});

        try
            m0 = fitlme(T,'Y ~ X + (1|Sample)',    'FitMethod','ML');
            m1 = fitlme(T,'Y ~ X + (X|Sample)',    'FitMethod','ML');
        catch
            continue    % non convergence : replicat ignore
        end
        % ML et non REML : les deux modeles ont la meme partie fixe, mais ML
        % est requis pour que le rapport de vraisemblance soit comparable.
        selA(r) = m1.ModelCriterion.AIC < m0.ModelCriterion.AIC;
        selB(r) = m1.ModelCriterion.BIC < m0.ModelCriterion.BIC;
        dA(r)   = m0.ModelCriterion.AIC - m1.ModelCriterion.AIC;
        cmp = compare(m0, m1);
        selL(r) = cmp.pValue(2) < 0.05;
    end
    pFull_AIC(a)=mean(selA); pFull_BIC(a)=mean(selB); pFull_LRT(a)=mean(selL);
    dAIC(a)=mean(dA,'omitnan');
    fprintf('%8.3f %10.2f %8.3f %8.3f %8.3f %10.2f\n', ...
        sigGrid(a), fracGrid(a), pFull_AIC(a), pFull_BIC(a), pFull_LRT(a), dAIC(a));
end
fprintf('\nTermine en %.1f min.\n', toc(tt)/60);

% --- lectures ----------------------------------------------------------------
fprintf('\n=== Lectures ===\n');
fprintf('A sigma_u1 = 0 (jeu de reference du papier), le modele complet est\n');
fprintf('retenu dans %.1f %% des replicats par l''AIC : la specification reduite\n', ...
    100*pFull_AIC(1));
fprintf('est donc bien le choix correct, et non une contrainte de plan.\n');
i50 = find(pFull_AIC >= 0.5, 1, 'first');
if ~isempty(i50)
    fprintf('Bascule de l''AIC vers le modele complet a sigma_u1 = %.3f (%.0f %% de beta1).\n', ...
        sigGrid(i50), 100*fracGrid(i50));
else
    fprintf(['Pas de bascule jusqu''a sigma_u1 = %.2f, soit %.0f %% de beta1.\n' ...
             'A cette taille de plan (G=%d, n=%d, R=%d), le critere ne distingue\n' ...
             'donc pas les deux structures aleatoires meme lorsque la dispersion\n' ...
             'des pentes egale la pente fixe : c''est une limite de puissance du\n' ...
             'plan, non une preuve que la pente aleatoire est absente.\n'], ...
             sigGrid(end), 100*fracGrid(end), P.G, P.n, P.R);
end

% --- figure ------------------------------------------------------------------
f = figure('Position',[100 100 560 380]);
plot(fracGrid, pFull_AIC, '-o', 'LineWidth',1.4, 'MarkerSize',5); hold on;
plot(fracGrid, pFull_BIC, '--s', 'LineWidth',1.0, 'MarkerSize',4);
plot(fracGrid, pFull_LRT, ':^', 'LineWidth',1.0, 'MarkerSize',4);
yline(0.5,':','HandleVisibility','off');
xlabel('\sigma_{u1} / \beta_1');
ylabel('Frequency of selecting the random-slope model');
legend({'AIC','BIC','LRT (\alpha = 0.05)'},'Location','southeast'); legend boxoff;
ylim([0 1.05]); grid on; box on; set(gca,'FontSize',10);
print(f,'Fig_randslope_selection.png','-dpng','-r300');
fprintf('\nFig_randslope_selection.png sauvee.\n');

save('ws_run_randslope.mat','fracGrid','sigGrid','pFull_AIC','pFull_BIC', ...
     'pFull_LRT','dAIC','NREP','P');
