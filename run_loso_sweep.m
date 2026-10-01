% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/08/2026 : calibration du drapeau fragile (section 4.7)
%RUN_LOSO_SWEEP  Balayage du facteur d'inflation de l'anomalie specimen.
%   run_fragile.m demontre que le drapeau fragile se leve, mais sur une seule
%   graine et une grille grossiere. Ce script etablit la courbe de calibration :
%   frequence de levee du drapeau en fonction du facteur d'inflation, sur
%   NREP replicats independants. C'est ce qui justifie la phrase du papier sur
%   l'amplitude minimale requise, aujourd'hui affirmee sans support.
%
%   Necessite : gen_fragile.m, gen_volet1.m, paper_pipeline.m, paper_params.m
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P = paper_params(); NREP = 200;

ampGrid = 1:20;                       % facteur d'inflation de la dispersion ROI
mech = 'variance'; sTgt = 5; kTgt = 3; iTgt = 1;

nA = numel(ampGrid);
freqAny = zeros(nA,1);    % au moins une paire fragile
freqTgt = zeros(nA,1);    % la paire cible elle-meme fragile
nFrag   = zeros(nA,1);    % nombre moyen de paires fragiles

fprintf('=== Calibration du drapeau fragile, %d replicats par point ===\n', NREP);
fprintf('mecanisme=%s, cible s=%d k=%d, specimen i=%d\n\n', mech, sTgt, kTgt, iTgt);

tt = tic;
for a = 1:nA
    amp = ampGrid(a);
    hitAny = false(NREP,1); hitTgt = false(NREP,1); cnt = zeros(NREP,1);
    parfor r = 1:NREP
        rs = RandStream('Threefry','Seed', paper_seed('fragile', r, 80+a));   % bloc 80-99
        [X,Y] = gen_fragile(P, rs, ...
            struct('mech',mech,'s',sTgt,'k',kTgt,'i',iTgt,'amp',amp));
        out = paper_pipeline(X,Y,P,struct('loso',true));
        hitAny(r) = nnz(out.fragile) >= 1;
        hitTgt(r) = out.fragile(sTgt,kTgt);
        cnt(r)    = nnz(out.fragile);
    end
    freqAny(a) = mean(hitAny); freqTgt(a) = mean(hitTgt); nFrag(a) = mean(cnt);
    fprintf('amp=%2d  drapeau (toute paire)=%.3f  (paire cible)=%.3f  nb moyen=%.2f\n', ...
        amp, freqAny(a), freqTgt(a), nFrag(a));
end
fprintf('\nTermine en %.1f min.\n', toc(tt)/60);

% --- seuil de declenchement --------------------------------------------------
i50 = find(freqTgt >= 0.50, 1, 'first');
i05 = find(freqTgt >= 0.05, 1, 'first');
if ~isempty(i05)
    fprintf('\nPremiere amplitude levant le drapeau dans >5%% des replicats  : %d\n', ampGrid(i05));
end
if ~isempty(i50)
    fprintf('Premiere amplitude levant le drapeau dans >50%% des replicats : %d\n', ampGrid(i50));
end

% --- figure ------------------------------------------------------------------
f = figure('Position',[100 100 560 380]);
plot(ampGrid, freqTgt, '-o', 'LineWidth',1.4, 'MarkerSize',5); hold on;
plot(ampGrid, freqAny, '--s', 'LineWidth',1.0, 'MarkerSize',4);
yline(0.5,':','HandleVisibility','off');
xlabel('ROI dispersion inflation factor of the anomalous specimen');
ylabel('Frequency of the fragility flag');
legend({'target pair','any pair'},'Location','southeast'); legend boxoff;
ylim([0 1.05]); xlim([min(ampGrid) max(ampGrid)]); grid on; box on;
set(gca,'FontSize',10);
print(f,'Fig_loso_sweep.png','-dpng','-r300');
fprintf('\nFig_loso_sweep.png sauvee.\n');

save('ws_run_loso_sweep.mat','ampGrid','freqAny','freqTgt','nFrag','NREP','P');
