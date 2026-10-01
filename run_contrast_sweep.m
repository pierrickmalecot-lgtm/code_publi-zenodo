% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 04/08/2026 : contraste entre conditions et grille de fiabilite
%RUN_CONTRAST_SWEEP  Ce qui gouverne la grille de fiabilite.
%
%   Le balayage de plan de la section 5.5 fait varier G a etendue de conditions
%   fixe. L'ICC y decroit avec G, ce qui se lit a tort comme une recommandation
%   d'utiliser peu de conditions. La cause est ailleurs : a etendue fixe,
%   ajouter des niveaux les rapproche, la variance entre conditions diminue, et
%   l'ICC suit. Ce n'est pas un defaut du plan, c'est le prix d'une question
%   plus fine.
%
%   Ce script separe les deux effets.
%     Volet A : ICC contre le contraste D/sigma, a G fixe. Abaque directement
%               utilisable : quel ecart entre la condition la plus lisse et la
%               plus rugueuse, rapporte a la dispersion intra-condition, est
%               necessaire pour franchir le seuil de 0.6.
%     Volet B : ICC contre G, a variance entre conditions egalisee. Si l'ICC y
%               est plat, le nombre de conditions n'a en lui-meme aucun effet
%               sur la fiabilite, et le resultat de la section 5.5 est bien un
%               artefact de parametrage.
%
%   D est l'etendue des niveaux de condition (le plus rugueux moins le plus
%   lisse) et sigma la dispersion totale intra-condition,
%   sigma = sqrt(sigma_between^2 + sigma_within^2), soit ce qu'un experimentateur
%   mesure sur des specimens nominalement identiques.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P = paper_params(); NREP = 200;

sigW   = sqrt(P.sigma_between^2 + P.sigma_within^2);   % dispersion intra-condition
s_star = P.s_star; kRef = 1; bandRef = 5:7;

fprintf('=== Contraste et fiabilite, %d replicats par point ===\n', NREP);
fprintf('dispersion intra-condition sigma = %.3f\n\n', sigW);

%% ---------------------------------------------------------------- volet A --
% ICC contre le contraste, a G = 5 fixe.
ratioGrid = [1 1.5 2 2.5 3 3.5 4 5 6 8];      % D / sigma
Dgrid     = ratioGrid * sigW;
nA = numel(Dgrid);
iccA = nan(nA,1); powA = nan(nA,1); r2A = nan(nA,1); ratioEff = nan(nA,1);

% Le contraste utile n'est pas celui des reglages du procede mais celui du
% DESCRIPTEUR a l'echelle d'analyse. Le generateur amplifie l'effet condition
% sur la bande (facteur 1 + 0.3(s - 5), soit 1.3 a s = 6), et sur une surface
% reelle un descripteur repond de meme avec sa propre sensibilite. Le rapport
% reporte ci-dessous est donc mesure sur X, comme le ferait un experimentateur
% sur une pre-campagne : ecart entre conditions extremes divise par l'ecart-type
% intra-condition.
fprintf('--- Volet A : contraste, G = %d ---\n', P.G);
fprintf('%8s %8s %8s %8s %8s\n','nominal','on X','ICC','Power','R2');
for a = 1:nA
    Pa = P;
    Pa.condLevels = linspace(0.2, 0.2+Dgrid(a), Pa.G);
    acc = nan(NREP,3); accRatio = nan(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('nominal', r, 20+a));   % bloc 20-39
        [X,Y] = gen_volet1(Pa,'nominal',rs);
        out   = paper_pipeline(X,Y,Pa);
        g1 = (out.p_fdr < Pa.q) & (out.eta2 >= Pa.eta2_thr);
        acc(r,:) = [out.ICC(s_star,kRef), mean(g1(bandRef,kRef)), ...
                    out.R2m(s_star,kRef)];
        % contraste effectif mesure sur le descripteur, a l'echelle d'analyse
        Xs = X(:,:,s_star,kRef);
        cm = arrayfun(@(g) mean(Xs(sampleCond(Pa)==g,:),'all'), 1:Pa.G);
        sw = sqrt(mean(arrayfun(@(g) var(reshape(Xs(sampleCond(Pa)==g,:),[],1),1), 1:Pa.G)));
        accRatio(r) = (max(cm)-min(cm))/sw;
    end
    m = mean(acc,1,'omitnan');
    iccA(a)=m(1); powA(a)=m(2); r2A(a)=m(3);
    ratioEff(a) = mean(accRatio,'omitnan');
    fprintf('%8.2f %8.2f %8.3f %8.3f %8.3f\n', ...
        ratioGrid(a), ratioEff(a), m(1), m(2), m(3));
end

% contraste minimal franchissant le seuil
iOK = find(iccA >= P.icc_thr, 1, 'first');
if ~isempty(iOK)
    if iOK > 1   % interpolation lineaire entre les deux points encadrants
        x1=ratioEff(iOK-1); x2=ratioEff(iOK); y1=iccA(iOK-1); y2=iccA(iOK);
        rStar = x1 + (P.icc_thr-y1)*(x2-x1)/(y2-y1);
    else
        rStar = ratioEff(1);
    end
    fprintf('\nContraste minimal sur le descripteur pour ICC >= %.2f : %.2f\n', P.icc_thr, rStar);
else
    fprintf('\nSeuil non atteint sur la grille : etendre ratioGrid.\n');
end

%% ---------------------------------------------------------------- volet B --
% ICC contre G, a variance entre conditions egalisee sur la valeur de G = 5.
Ggrid = [3 5 7 9];
vRef  = var(linspace(0, 1, 5), 1) * (0.8)^2;   % variance des niveaux de reference
nB = numel(Ggrid);
iccB = nan(nB,1); powB = nan(nB,1); fpB = nan(nB,1); Dused = nan(nB,1);

fprintf('\n--- Volet B : nombre de conditions, variance entre conditions egalisee ---\n');
fprintf('%4s %8s %8s %8s %8s\n','G','D','ICC','Power','TypeI');
for b = 1:nB
    G = Ggrid(b);
    k = var(linspace(0,1,G), 1);          % var pour une etendue unitaire
    D = sqrt(vRef / k);                   % etendue donnant la variance de reference
    Dused(b) = D;
    Pb = P; Pb.G = G; Pb.N = G*Pb.n;
    Pb.condLevels = linspace(0.2, 0.2+D, G);
    acc = nan(NREP,3);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('nominal', r, 40+b));   % bloc 40-59
        [X,Y] = gen_volet1(Pb,'nominal',rs);
        out   = paper_pipeline(X,Y,Pb);
        g1 = (out.p_fdr < Pb.q) & (out.eta2 >= Pb.eta2_thr);
        acc(r,:) = [out.ICC(s_star,kRef), mean(g1(bandRef,kRef)), mean(g1(:,4))];
    end
    m = mean(acc,1,'omitnan');
    iccB(b)=m(1); powB(b)=m(2); fpB(b)=m(3);
    fprintf('%4d %8.3f %8.3f %8.3f %8.4f\n', G, D, m(1), m(2), m(3));
end
fprintf('\nEcart max d''ICC sur le volet B : %.3f\n', max(iccB)-min(iccB));
fprintf('Si cet ecart est faible, le nombre de conditions n''a pas d''effet propre\n');
fprintf('sur la fiabilite, et la decroissance de la section 5.5 est un artefact.\n');

%% ------------------------------------------------------------------ figure --
f = figure('Position',[100 100 900 360]);

subplot(1,2,1); hold on;
plot(ratioEff, iccA, '-o','LineWidth',1.5,'MarkerSize',5);
% Courbe de reference ICC = c r^2 / (c r^2 + 1), la constante c etant estimee
% par moindres carres sur les points simules plutot que derivee des niveaux
% nominaux : le descripteur repond a la condition avec sa propre sensibilite,
% inconnue a priori sur une surface reelle.
cHat = mean( iccA(:)./(1-iccA(:)) ./ (ratioEff(:).^2) );
rr = linspace(min(ratioEff), max(ratioEff), 200);
plot(rr, cHat*(rr.^2)./(cHat*(rr.^2)+1), '--','LineWidth',1.0);
yline(P.icc_thr,':','HandleVisibility','off');
xlabel('Descriptor contrast  \Delta / \sigma_{within}');
ylabel('ICC'); ylim([0 1]); grid on; box on;
legend({'simulated','fitted'},'Location','southeast'); legend boxoff;
title('Reliability against contrast');

subplot(1,2,2); hold on;
plot(Ggrid, iccB, '-o','LineWidth',1.5,'MarkerSize',5);
plot(Ggrid, powB, '-s','LineWidth',1.5,'MarkerSize',5);
yline(P.icc_thr,':','HandleVisibility','off');
xlabel('Number of conditions G'); ylabel('Value'); ylim([0 1.05]);
xticks(Ggrid); grid on; box on;
legend({'ICC','Power'},'Location','southeast'); legend boxoff;
title('At matched between-condition variance');

print(f,'Fig_contrast_sweep.png','-dpng','-r300');
fprintf('\nFig_contrast_sweep.png sauvee.\n');

save('ws_run_contrast_sweep.mat','ratioGrid','ratioEff','Dgrid','iccA','powA','r2A', ...
     'Ggrid','Dused','iccB','powB','fpB','sigW','NREP','P');


function g = sampleCond(P)
%SAMPLECOND  Condition de chaque specimen, meme convention que gen_volet1.
g = repelem((1:P.G)', P.n, 1);
end
