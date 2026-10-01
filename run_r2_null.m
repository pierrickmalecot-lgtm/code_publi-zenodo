% version v1, 29/09/2026 : R2 marginal obtenu par hasard, sans lien X-Y
%RUN_R2_NULL  Distribution du R2 marginal quand la reponse ne depend pas du
%   descripteur (beta1 = 0), en fonction du nombre de specimens par condition.
%
%   Deux lectures par replicat :
%     - la paire porteuse D1 a s* (une paire fixee a l'avance) ;
%     - le maximum sur les S x K paires, c'est-a-dire la meilleure paire que
%       trouverait un criblage par correlation seule, sans verrous.
%   Generateur 'nominal', pente vraie mise a zero, R = 5 fixe.
%
%   Graines : bloc 510-519 (a declarer dans paper_seed).
%   Necessite : paper_params, gen_volet1, paper_seed.

P0    = paper_params();
NREP  = 200;
ngrid = [3 5 10];
nn    = numel(ngrid);
S = P0.S; K = P0.K; s = P0.s_star; k = 1;

R2pair = nan(NREP,nn);      % D1 a s*
R2max  = nan(NREP,nn);      % max sur toutes les paires
for in = 1:nn
    P = P0; P.n = ngrid(in); P.N = P.G*P.n; P.beta1 = 0;
    lev = 509 + in;
    rp = nan(NREP,1); rm = nan(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('nominal', r, lev));
        [X,Y] = gen_volet1(P,'nominal',rs);
        v = nan(S,K);
        for ss = 1:S
            for kk = 1:K
                try, v(ss,kk) = r2_marginal(X(:,:,ss,kk), Y); catch, end
            end
        end
        rp(r) = v(s,k);
        rm(r) = max(v(:));
    end
    R2pair(:,in) = rp; R2max(:,in) = rm;
end

% --- synthese ---------------------------------------------------------------
fprintf('\nR2 marginal sans lien X-Y, %d replicats, R = %d\n', NREP, P0.R);
fprintf('%4s | %28s | %28s\n','n','D1 a s* : med / p95 / max','max sur paires : med / p95 / max');
res = nan(nn,7);
for in = 1:nn
    a = R2pair(:,in); b = R2max(:,in);
    q = [median(a,'omitnan') prctile(a,95) max(a) median(b,'omitnan') prctile(b,95) max(b)];
    fprintf('%4d | %8.3f %8.3f %8.3f   | %8.3f %8.3f %8.3f\n', ngrid(in), q);
    res(in,:) = [ngrid(in) q];
end
fprintf('Part des replicats ou le max sur paires depasse %.2f : ', P0.r2_thr);
fprintf('%.3f  ', mean(R2max > P0.r2_thr, 1)); fprintf('\n');
writematrix(res, 'r2_null.csv');
save('ws_run_r2_null.mat','R2pair','R2max','ngrid','res');

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
