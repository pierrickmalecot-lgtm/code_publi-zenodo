% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/08/2026 : courbure x bruit, separation des deux mecanismes
%RUN_CURVATURE_NOISE  Balayage croise du taux de courbure et du niveau de bruit
%   de la reponse.
%
%   Objet : etablir que les deux indices de correlation se degradent pour des
%   raisons distinctes. Le R2 marginal souffre de la courbure, qui envoie une
%   part de la relation dans le residu, et cette perte ne devrait pas dependre
%   du bruit. Spearman ne souffre pas de la courbure elle-meme, la
%   transformation restant strictement monotone, mais du bruit dans la zone ou
%   la reponse sature et ou l'ordre des specimens n'est plus gouverne par la
%   topographie ; sa perte devrait donc s'attenuer quand le bruit diminue.
%
%   Lecture attendue : le long d'une colonne (bruit fixe) les deux indices
%   baissent ; le long d'une ligne (courbure fixe) seul rho remonte quand le
%   bruit baisse. La perte est mesuree relativement a la reference lineaire du
%   MEME niveau de bruit, faute de quoi la comparaison serait faussee par
%   l'augmentation mecanique du R2 a bruit reduit.
%
%   Necessite gen_volet1.m v7 (scenario 'curvature') et paper_params.m v7.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P = paper_params(); NREP = 200;

curvGrid  = [0 1 3 5];                 % 0 = reference lineaire
noiseGrid = [0.25 0.50 1.00];          % facteur applique a sigma_func et sigma_u
noiseLbl  = {'low','half','reference'};

s_star = P.s_star; kRef = 1;
nC = numel(curvGrid); nN = numel(noiseGrid);
R2  = nan(nC,nN); RHO = nan(nC,nN); ICC = nan(nC,nN);

fprintf('=== Courbure x bruit, %d replicats par case ===\n', NREP);
fprintf('sigma_func = %.2f, sigma_u = %.2f a facteur 1\n\n', P.sigma_func, P.sigma_u);

tt=tic;
for j = 1:nN
    for i = 1:nC
        Pc = P;
        Pc.curv       = curvGrid(i);
        Pc.sigma_func = P.sigma_func * noiseGrid(j);
        Pc.sigma_u    = P.sigma_u    * noiseGrid(j);

        acc = nan(NREP,3);
        parfor mc = 1:NREP
            rs = RandStream('Threefry','Seed', paper_seed('curvature', mc, 140+10*j+i));   % bloc 140-199
            [X,Y] = gen_volet1(Pc,'curvature',rs);
            out   = paper_pipeline(X,Y,Pc);
            acc(mc,:) = [out.R2m(s_star,kRef), out.rho(s_star,kRef), ...
                         out.ICC(s_star,kRef)];
        end
        m = mean(acc,1,'omitnan');
        R2(i,j)=m(1); RHO(i,j)=m(2); ICC(i,j)=m(3);
        fprintf('bruit %-9s c=%d : R2=%.3f  rho=%.3f  ICC=%.3f\n', ...
            noiseLbl{j}, curvGrid(i), m(1), m(2), m(3));
    end
    fprintf('\n');
end
fprintf('Termine en %.1f min.\n', toc(tt)/60);

% --- pertes relatives, chaque colonne rapportee a sa propre reference --------
lossR2  = 100*(R2  - R2(1,:) ) ./ R2(1,:);
lossRHO = 100*(RHO - RHO(1,:)) ./ RHO(1,:);

fprintf('\n=== Perte relative (%%) par rapport a la reference lineaire du meme bruit ===\n');
fprintf('%-10s', 'courbure');
for j=1:nN, fprintf('%18s', noiseLbl{j}); end
fprintf('\n');
for i=2:nC
    fprintf('c = %-6d', curvGrid(i));
    for j=1:nN
        fprintf('   R2 %+6.1f rho %+6.1f', lossR2(i,j), lossRHO(i,j));
    end
    fprintf('\n');
end

fprintf('\nSi la these tient : la perte de R2 varie peu le long d''une ligne,\n');
fprintf('celle de rho s''attenue nettement quand le bruit diminue.\n');

% --- figure ------------------------------------------------------------------
f = figure('Position',[100 100 900 360]);
mk = {'-o','-s','-^'};
subplot(1,2,1); hold on;
for j=1:nN, plot(curvGrid, R2(:,j), mk{j}, 'LineWidth',1.3,'MarkerSize',5); end
xlabel('Curvature rate c'); ylabel('R^2_{marginal}');
legend(noiseLbl,'Location','southwest'); legend boxoff;
grid on; box on; ylim([0 1]); title('Linear fit');
subplot(1,2,2); hold on;
for j=1:nN, plot(curvGrid, RHO(:,j), mk{j}, 'LineWidth',1.3,'MarkerSize',5); end
xlabel('Curvature rate c'); ylabel('Spearman \rho');
grid on; box on; ylim([0 1]); title('Rank association');
print(f,'Fig_curvature_noise.png','-dpng','-r300');
fprintf('\nFig_curvature_noise.png sauvee.\n');

save('ws_run_curvature_noise.mat','curvGrid','noiseGrid','noiseLbl', ...
     'R2','RHO','ICC','lossR2','lossRHO','NREP','P');
