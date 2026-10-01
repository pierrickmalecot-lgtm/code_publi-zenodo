% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/09/2026 : effet du nombre de specimens sur les tests permutes
%RUN_PERM_SIZE  A partir de combien de specimens la permutation studentisee
%   tient-elle le niveau nominal sous variances inegales.
%
%   CE QUE CE SCRIPT TESTE. La mesure precedente donne, sous heteroscedasticite
%   et sur D4, un taux d'erreur de type I de 0.0735 pour la permutation
%   studentisee de Welch, soit le meme que pour la permutation du F ordinaire
%   et au-dessus du 0.0615 du logigramme. Or studentiser la statistique est
%   cense rendre le test valide sous variances inegales. L'explication
%   proposee est une question de taille : la validite des tests permutes
%   studentises est un resultat ASYMPTOTIQUE en la taille des groupes, et
%   l'unite de permutation est ici le SPECIMEN, soit n = 5 unites par condition
%   et 25 en tout. C'est tres loin du regime asymptotique.
%
%   Si cette explication est juste, l'ecart au niveau nominal doit decroitre
%   quand n augmente, pour la permutation de Welch et pour elle seule. Si
%   l'ecart persiste a n = 40, l'explication est fausse et il faut en chercher
%   une autre.
%
%   ATTENTION A CE QUI VARIE. Seul n change. G, R, S, K et toutes les
%   amplitudes de bruit restent aux valeurs du papier, donc les rapports de
%   variance entre conditions du scenario heteroscedastique sont inchanges.
%   Le nombre d'observations poolees par condition, m = n*R, augmente avec n :
%   c'est voulu, c'est exactement la quantite dont depend la statistique de
%   Welch.
%
%   DEUX NULLS SUIVIS. D4, bruit de region pur, mesure le niveau sous variances
%   inegales. D5, bruit specimen plus region, verifie au passage que la
%   correction de la pseudo-replication par la permutation specimen tient a
%   toutes les tailles.
%
%   COUT. Il croit vite avec n : le verrou 3 ajuste des modeles mixtes sur
%   n*G*R lignes. NREP est donc reduit pour les grandes valeurs de n, au prix
%   d'une incertitude plus large, indiquee dans la sortie.
%
%   Necessite : gen_volet1, paper_pipeline v8, paper_params.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NGRID  = [5 10 20 40];        % specimens par condition
NREP_N = [200 200 100 100];   % replicats, reduits pour les grandes tailles
B_PERM = 10000;
SCEN   = 'heteroscedastic';
METH   = {'branching','perm','permwelch'};
MLBL   = {'logigramme papier','permutation (F)','permutation (Welch)'};
KNULL  = [4 5];

P0 = paper_params();
S = P0.S; NM = numel(METH); NN = numel(NGRID); NK = numel(KNULL);

alpha_null = zeros(NN, NM, NK);
fp_null    = zeros(NN, NM, NK);
band_D1    = zeros(NN, NM);     % paires D1 au verrou 1 ; la bande vraie fait 3 echelles
n_gate1    = zeros(NN, NM);
n_gate3    = zeros(NN, NM);
se_alpha   = zeros(NN, 1);      % incertitude sur alpha, ecart-type binomial

fprintf('=== Effet du nombre de specimens, scenario %s, B_perm = %d ===\n\n', ...
    SCEN, B_PERM);
tt = tic;

for in = 1:NN
    nn = NGRID(in); NREP = NREP_N(in);
    P = P0; P.n = nn; P.N = P.G*nn;      % seule la taille change
    se_alpha(in) = sqrt(0.05*0.95/(S*NREP));

    aN = zeros(NREP,NM,NK); fN = aN;
    bD = zeros(NREP,NM); g1 = bD; g3 = bD;

    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed(SCEN, r, 100+in));   % bloc 100-119
        [X,Y] = gen_volet1(P,SCEN,rs);
        aN_r = zeros(NM,NK); fN_r = aN_r;
        bD_r = zeros(1,NM); g1_r = bD_r; g3_r = bD_r;
        for im = 1:NM
            o = paper_pipeline(X,Y,P, struct('gate1Test',METH{im}, ...
                                             'B_perm',B_PERM));
            gate1 = (o.p_fdr < P.q) & (o.eta2 >= P.eta2_thr);
            for ik = 1:NK
                kk = KNULL(ik);
                aN_r(im,ik) = mean(o.p_raw(:,kk) < 0.05);
                fN_r(im,ik) = sum(gate1(:,kk));
            end
            bD_r(im) = sum(gate1(:,1));
            g1_r(im) = sum(gate1(:));
            g3_r(im) = sum(o.decision(:) == 3);
        end
        aN(r,:,:)=aN_r; fN(r,:,:)=fN_r;
        bD(r,:)=bD_r; g1(r,:)=g1_r; g3(r,:)=g3_r;
    end

    alpha_null(in,:,:) = squeeze(mean(aN,1));
    fp_null(in,:,:)    = squeeze(mean(fN,1));
    band_D1(in,:)      = mean(bD,1);
    n_gate1(in,:)      = mean(g1,1);
    n_gate3(in,:)      = mean(g3,1);
    fprintf('n = %2d (%3d replicats) termine, %.1f min ecoulees\n', ...
        nn, NREP, toc(tt)/60);
end
fprintf('\nTotal %.1f min.\n', toc(tt)/60);

% --- sorties ------------------------------------------------------------------
fprintf('\n--- D4, bruit de region pur : niveau sous variances inegales ---\n');
fprintf('Niveau nominal 0.05. L incertitude sur alpha est donnee par taille.\n');
fprintf('%4s %8s %6s   %-22s  %10s  %14s\n', ...
    'n', 'm=n*R', '+/-', 'methode', 'alpha_D4', 'paires verrou 1');
for in = 1:NN
    for im = 1:NM
        fprintf('%4d %8d %6.3f   %-22s  %10.4f  %14.3f\n', ...
            NGRID(in), NGRID(in)*P0.R, se_alpha(in), MLBL{im}, ...
            alpha_null(in,im,1), fp_null(in,im,1));
    end
end

fprintf('\n--- D5, bruit specimen + region : pseudo-replication ---\n');
fprintf('%4s   %-22s  %10s  %14s\n','n','methode','alpha_D5','paires verrou 1');
for in = 1:NN
    for im = 1:NM
        fprintf('%4d   %-22s  %10.4f  %14.3f\n', NGRID(in), MLBL{im}, ...
            alpha_null(in,im,2), fp_null(in,im,2));
    end
end

fprintf('\n--- Localisation de la bande de D1 (vraie largeur : 3 echelles sur %d) ---\n', S);
fprintf('%4s   %-22s  %10s  %12s  %12s\n', ...
    'n','methode','D1 verrou 1','total verrou 1','retenues');
for in = 1:NN
    for im = 1:NM
        fprintf('%4d   %-22s  %10.2f  %12.2f  %12.2f\n', NGRID(in), MLBL{im}, ...
            band_D1(in,im), n_gate1(in,im), n_gate3(in,im));
    end
end

% --- figure -------------------------------------------------------------------
% Les libelles de la legende sont en anglais, ceux de MLBL servant a l'affichage
% console en francais. La legende est placee explicitement et non en 'best', qui
% la posait sur une courbe.
LEG = {'branching','permutation (F)','permutation (Welch)'};
f = figure;
mk = {'-o','-s','-^'};
for panel = 1:2
    subplot(1,2,panel); hold on;
    for im = 1:NM
        plot(NGRID, alpha_null(:,im,panel), mk{im}, 'LineWidth',1.2, 'MarkerSize',5);
    end
    yline(0.05,'--','HandleVisibility','off');
    set(gca,'XScale','log','XTick',NGRID,'XTickLabel',string(NGRID));
    xlabel('Specimens per condition, n');
    ylabel('Type-I error rate');
    if panel==1
        title('D4: pure region noise');
        legend(LEG,'Location','northeast');
    else
        title('D5: specimen + region noise');
    end
    grid on; box on;
end
paper_figstyle(f, 15, 7);
print(f,'Fig_perm_size.png','-dpng','-r300');
fprintf('\nFig_perm_size.png sauvee.\n');

save('ws_run_perm_size.mat','alpha_null','fp_null','band_D1','n_gate1', ...
     'n_gate3','se_alpha','NGRID','NREP_N','METH','MLBL','KNULL','B_PERM','SCEN','P0');
