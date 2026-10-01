% version v2, 05/09/2026 : texte de lecture corrige, le resultat est l inverse
%                          de l hypothese initiale
% version v1, 04/09/2026 : ce qui restaure la puissance sur la pente aleatoire
%RUN_RANDSLOPE_R  Le nombre de regions par specimen, ou le nombre de specimens.
%
%   CE QUI PRECEDE. run_randslope balaie sigma_u1 de 0 a 1.0*beta1 sur le plan
%   de reference, G=5, n=5, R=5. Le critere d'information ne bascule jamais
%   vers le modele a pente aleatoire : la proportion de replicats ou l'AIC le
%   retient passe de 4 % a environ 17 %, sans tendance nette au-dela de 0.2, et
%   l'ecart d'AIC reste negatif partout. La specification reduite n'est donc
%   pas SELECTIONNEE par les donnees, elle est conservee faute que le critere
%   puisse trancher.
%
%   HYPOTHESE TESTEE ICI. La pente d'un specimen est estimee sur ses R regions.
%   L'information sur la dispersion des pentes est portee par la variation
%   INTRA-specimen, pas par le nombre de specimens : ajouter des specimens
%   ajoute des pentes toutes aussi mal estimees les unes que les autres, tandis
%   qu'ajouter des regions estime mieux chacune. Si l'hypothese est juste, un
%   balayage en R doit restaurer la puissance et un balayage en n beaucoup
%   moins, a nombre total d'observations comparable.
%
%   PROTOCOLE. sigma_u1 est fixe a une valeur franche, 0.6*beta1, ou la pente
%   aleatoire existe sans ambiguite. Deux balayages, chacun a l'autre facteur
%   fixe a sa valeur de reference :
%     volet A : R = 5, 10, 20, 40 a n = 5
%     volet B : n = 5, 10, 20, 40 a R = 5
%   Les deux atteignent le meme nombre total d'observations aux memes rangs,
%   ce qui rend la comparaison lisible : si seul le volet A monte, l'effet
%   vient bien des regions et non de la taille de l'echantillon.
%
%   Necessite : gen_volet1 (scenario 'randslope'), paper_params, paper_seed.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP = 200;
FRAC = 0.60;          % sigma_u1 = FRAC * beta1, pente aleatoire franche
GRID = [5 10 20 40];  % memes valeurs pour les deux volets

P0 = paper_params();
s_star = P0.s_star; kRef = 1;
nG = numel(GRID);

res = struct();
for volet = 1:2
    if volet==1, nom = 'A : regions par specimen (R)';
    else,        nom = 'B : specimens par condition (n)'; end
    pA = zeros(nG,1); pB = pA; pL = pA; dA = pA; nobs = pA;

    for a = 1:nG
        Pa = P0; Pa.sigma_u1 = FRAC * P0.beta1;
        if volet==1
            Pa.R = GRID(a);                       % n reste a sa valeur de reference
        else
            Pa.n = GRID(a); Pa.N = Pa.G*Pa.n;     % R reste a sa valeur de reference
        end
        nobs(a) = Pa.N * Pa.R;
        lvl = 220 + 10*volet + a;                 % blocs 230-239 et 240-249

        selA = false(NREP,1); selB = selA; selL = selA; d = nan(NREP,1);
        parfor r = 1:NREP
            rs = RandStream('Threefry','Seed', paper_seed('randslope', r, lvl));
            [X,Y] = gen_volet1(Pa,'randslope',rs);
            x = reshape(X(:,:,s_star,kRef),[],1); x = (x-mean(x))/std(x);
            y = Y(:);
            sid = repelem((1:Pa.N)', Pa.R, 1);
            T = table(y, x, nominal(sid), 'VariableNames',{'Y','X','Sample'});
            try
                m0 = fitlme(T,'Y ~ X + (1|Sample)', 'FitMethod','ML');
                m1 = fitlme(T,'Y ~ X + (X|Sample)', 'FitMethod','ML');
            catch
                continue
            end
            selA(r) = m1.ModelCriterion.AIC < m0.ModelCriterion.AIC;
            selB(r) = m1.ModelCriterion.BIC < m0.ModelCriterion.BIC;
            d(r)    = m0.ModelCriterion.AIC - m1.ModelCriterion.AIC;
            cmp = compare(m0, m1);
            selL(r) = cmp.pValue(2) < 0.05;
        end
        pA(a)=mean(selA); pB(a)=mean(selB); pL(a)=mean(selL);
        dA(a)=mean(d,'omitnan');
    end
    res(volet).nom=nom; res(volet).pA=pA; res(volet).pB=pB;
    res(volet).pL=pL; res(volet).dA=dA; res(volet).nobs=nobs;
end

% --- sorties ------------------------------------------------------------------
fprintf('=== Pente aleatoire : ce qui restaure la puissance ===\n');
fprintf('sigma_u1 = %.2f (%.0f %% de beta1), %d replicats par point\n\n', ...
    FRAC*P0.beta1, 100*FRAC, NREP);
for volet = 1:2
    fprintf('--- Volet %s ---\n', res(volet).nom);
    fprintf('%8s %10s %8s %8s %8s %10s\n', ...
        'valeur','N x R','AIC','BIC','LRT','dAIC');
    for a = 1:nG
        fprintf('%8d %10d %8.3f %8.3f %8.3f %10.2f\n', GRID(a), ...
            res(volet).nobs(a), res(volet).pA(a), res(volet).pB(a), ...
            res(volet).pL(a), res(volet).dA(a));
    end
    fprintf('\n');
end

fprintf(['LECTURE. La colonne AIC donne la proportion de replicats ou le\n' ...
         'modele a pente aleatoire est retenu. Il devrait l etre presque\n' ...
         'toujours, la pente aleatoire etant reellement presente.\n' ...
         'Le resultat mesure est que le volet B monte et le volet A non :\n' ...
         'ce que le modele doit estimer est la VARIANCE des pentes entre\n' ...
         'specimens, un parametre de population dont la precision depend du\n' ...
         'nombre de specimens tires, et non de la precision avec laquelle\n' ...
         'chaque pente individuelle est connue. Ajouter des regions resserre\n' ...
         'chaque pente sans rien dire de leur dispersion : au volet A l ecart\n' ...
         'd AIC converge vers -4.00, soit exactement la penalite des deux\n' ...
         'parametres supplementaires, le modele etendu n achetant aucune\n' ...
         'vraisemblance. Un plan visant a tester cette structure doit donc\n' ...
         'multiplier les specimens, non les mesures faites sur chacun.\n']);

% --- figure -------------------------------------------------------------------
f = figure('Units','centimeters','Position',[2 2 15 7],'Color','w');
subplot(1,2,1); hold on;
plot(GRID, res(1).pA, '-o', 'LineWidth',1.4, 'MarkerSize',5);
plot(GRID, res(2).pA, '--s','LineWidth',1.4, 'MarkerSize',5);
yline(0.5,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',GRID,'XTickLabel',string(GRID),'FontSize',9);
xlabel('Regions per specimen, or specimens per condition','FontSize',9);
ylabel('P(AIC selects random slope)','FontSize',9);
legend({'sweep in R','sweep in n'},'Location','northwest','Box','off','FontSize',8);
grid on; box on; ylim([0 1]);

subplot(1,2,2); hold on;
plot(GRID, res(1).dA, '-o', 'LineWidth',1.4, 'MarkerSize',5);
plot(GRID, res(2).dA, '--s','LineWidth',1.4, 'MarkerSize',5);
yline(0,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',GRID,'XTickLabel',string(GRID),'FontSize',9);
xlabel('Regions per specimen, or specimens per condition','FontSize',9);
ylabel('\DeltaAIC (reduced - full)','FontSize',9);
grid on; box on;

print(f,'Fig_randslope_R.png','-dpng','-r300');
fprintf('\nFig_randslope_R.png sauvee.\n');

% La figure n'est pas sauvegardee dans le .mat : un objet graphique y occupe
% plusieurs mega-octets et MATLAB emet un avertissement a chaque enregistrement.
save('ws_run_randslope_R.mat','res','GRID','FRAC','NREP','P0');
