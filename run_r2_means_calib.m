% version v1, 28/09/2026 : calibration du seuil R2 du protocole simplifie
%RUN_R2_MEANS_CALIB  Seuil sur le R2 de l'OLS sur moyennes par specimen
%   equivalent au seuil P.r2_thr applique au R2 marginal du LME.
%
%   Paire porteuse D1 (k = 1) a l'echelle s*, generateur 'nominal'. La pente
%   vraie est balayee de 0 a P.beta1 pour couvrir tout l'intervalle de R2 :
%   au seul point nominal, les R2 sont tous loin du seuil et ne diraient rien.
%   Le balayage ne modifie que Y ; X, donc les verrous 1 et 2, sont inchanges.
%   Le calcul est repete pour R = 3, 5, 10, l'equivalence dependant de R.
%
%   Graines : bloc 460-504 (a declarer dans paper_seed).
%   Necessite : paper_params, gen_volet1, paper_seed.

P0    = paper_params();
NREP  = 100;
b1grid = linspace(0, P0.beta1, 15);
Rgrid  = [3 5 10];
k = 1; s = P0.s_star;
nb = numel(b1grid); nR = numel(Rgrid);

A = cell(nR,1); B = cell(nR,1);          % R2 marginal, R2 sur moyennes
for iR = 1:nR
    a = nan(NREP,nb); b = nan(NREP,nb);
    for ib = 1:nb
        P = P0; P.R = Rgrid(iR); P.beta1 = b1grid(ib);
        lev = 460 + (iR-1)*nb + ib - 1;
        ra = nan(NREP,1); rb = nan(NREP,1);
        parfor r = 1:NREP
            rs = RandStream('Threefry','Seed', paper_seed('nominal', r, lev));
            [X,Y] = gen_volet1(P,'nominal',rs);
            Xsk = X(:,:,s,k);
            try, ra(r) = r2_marginal(Xsk, Y); catch, end
            c = corrcoef(mean(Xsk,2), mean(Y,2));
            rb(r) = c(1,2)^2;
        end
        a(:,ib) = ra; b(:,ib) = rb;
    end
    A{iR} = a(:); B{iR} = b(:);
end

% --- seuil sur moyennes reproduisant les decisions du seuil marginal -------
tgrid = 0:0.01:1;
fprintf('\nSeuil marginal de reference : R2 > %.2f\n', P0.r2_thr);
fprintf('%4s %8s %8s %10s %14s %8s\n','R','t*','accord','plage t','R2moy a 0.3','n pres');
res = nan(nR,6);
for iR = 1:nR
    a = A{iR}; b = B{iR}; ok = ~isnan(a) & ~isnan(b); a = a(ok); b = b(ok);
    ref = a > P0.r2_thr;
    agr = arrayfun(@(t) mean(ref == (b > t)), tgrid);
    [amax, it] = max(agr);
    tOK = tgrid(agr >= amax - 0.005);
    near = abs(a - P0.r2_thr) < 0.03;           % paires proches du seuil
    fprintf('%4d %8.2f %8.3f %4.2f-%4.2f %14.3f %8d\n', Rgrid(iR), tgrid(it), amax, ...
        min(tOK), max(tOK), median(b(near)), nnz(near));
    res(iR,:) = [Rgrid(iR) tgrid(it) amax min(tOK) max(tOK) median(b(near))];
end
writematrix(res, 'r2_means_calib.csv');

% --- figure -----------------------------------------------------------------
f = figure('Position',[100 100 900 320]);
for iR = 1:nR
    subplot(1,nR,iR);
    plot(A{iR}, B{iR}, '.', 'MarkerSize', 4); hold on;
    plot([0 1],[0 1],'k:'); xline(P0.r2_thr,'r--'); yline(res(iR,2),'b--');
    axis([0 1 0 1]); axis square; grid on;
    xlabel('Marginal R^2 (LME)'); ylabel('R^2 on specimen means (OLS)');
    title(sprintf('R = %d, t^* = %.2f', Rgrid(iR), res(iR,2)));
end
print(f,'Fig_r2_means_calib.png','-dpng','-r300');
save('ws_run_r2_means_calib.mat','A','B','res','b1grid','Rgrid');

% ----------------------------------------------------------------------------
function R2m = r2_marginal(Xsk, Y)
% Copie de lme_r2_marginal (paper_pipeline) : memes conventions.
[N,R] = size(Xsk);
x = Xsk(:); x = (x-mean(x))/std(x); y = Y(:);
sid = repelem((1:N)',R,1);
T = table(y,x,nominal(sid),'VariableNames',{'Y','X','Sample'});
lme = fitlme(T,'Y ~ X + (1|Sample)','FitMethod','REML');
b1 = lme.Coefficients.Estimate(2);
varF = var(b1*x);
vc = covarianceParameters(lme); varU = vc{1}(1,1);
varE = lme.MSE;
R2m = varF/(varF+varU+varE);
end
