% version v9, 27/09/2026 : ecart-type bootstrap corrige n/(n-1) (o.b1_se du
%                         pipeline v10) pour le rapport ; biais normalise
%                         inchange par la correction
% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v1, 02/09/2026 : ce que le bootstrap capte et que le modele ignore
%RUN_BOOT_VS_PARAM  Ecart-type bootstrap de la pente contre ecart-type
%   parametrique, sur replicats.
%
%   POURQUOI. Sur le tirage de reference, la coherence de signe C vaut 1.000
%   sur les neuf paires retenues : elle sature des que l'intervalle ne traverse
%   pas zero, ce qui est le cas de toute paire ayant franchi les deux verrous
%   amont. C est donc constante sur son domaine d'emploi et ne hierarchise
%   rien. Dans le meme temps la precision relative varie d'un facteur quatre
%   entre paires retenues, et l'ecart-type bootstrap depasse l'ecart-type
%   parametrique d'un facteur 1.1 a 2.1 selon la paire.
%
%   CE QUE MESURE LE RAPPORT. L'ecart-type rendu par fitlme est calcule a
%   echantillon de specimens FIXE : il ne compte que le bruit residuel et la
%   dispersion des interceptes aleatoires. Le bootstrap stratifie reechantillonne
%   les specimens, il ajoute donc la variabilite due au fait que ces n specimens
%   par condition ne sont qu'un tirage parmi les possibles. Un rapport de 2
%   signifie que la moitie de l'incertitude sur la pente vient du choix des
%   specimens et non de la mesure : sur cette paire, ajouter des regions de
%   mesure ne servirait a rien, il faut plus de specimens.
%
%   CE QUE CE SCRIPT ETABLIT.
%     (1) Le rapport est-il stable sur replicats, ou l'ordre observe sur un
%         tirage tient-il a ce tirage.
%     (2) Le biais bootstrap, defini comme moyenne des pentes moins pente
%         estimee et rapporte a l'ecart-type bootstrap, est-il reproductible.
%         S'il l'est, il se corrige (intervalle de base ou BCa). S'il varie
%         d'un replicat a l'autre, aucune correction ne s'applique et
%         l'intervalle percentile doit etre abandonne sans remplacant.
%
%   COUT. Le bootstrap tourne sur TOUTES les paires retenues, soit environ neuf
%   paires par replicat. B_BOOT est donc reduit a 200 : ce qui est estime ici
%   est un ecart-type, pas un quantile extreme, et 200 tirages y suffisent.
%   Reduire NREP en premier si la duree est excessive.
%
%   Necessite : gen_volet1, paper_pipeline v4 (champs b1_se et Cboot),
%   paper_params.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('Processes'); end

NREP   = 200;      % replicats
B_BOOT = 200;      % tirages bootstrap par paire
FREQ_MIN = 0.20;   % frequence de retention minimale pour figurer dans la sortie

P = paper_params();
S = P.S; K = P.K;
names = {'D1','D2','D3','D4','D5'};

ratio = nan(NREP,S,K);   % se_boot / se_par
biasN = nan(NREP,S,K);   % (moyenne bootstrap - pente estimee) / se_boot
skw   = nan(NREP,S,K);   % asymetrie de la distribution bootstrap
seB   = nan(NREP,S,K);   % ecart-type bootstrap
seP   = nan(NREP,S,K);   % ecart-type parametrique
kept  = false(NREP,S,K); % paire retenue au verrou 3

fprintf('=== Bootstrap contre parametrique, %d replicats, B = %d ===\n\n', ...
    NREP, B_BOOT);
tt = tic;

parfor r = 1:NREP
    rs = RandStream('Threefry','Seed', paper_seed('nominal', r));
    [X,Y] = gen_volet1(P,'nominal',rs);
    o = paper_pipeline(X,Y,P,struct('bootstrap',true,'B',B_BOOT));

    rt = nan(S,K); bn = nan(S,K); sk = nan(S,K);
    sb = nan(S,K); sp = nan(S,K); kp = false(S,K);

    [si,ki] = find(o.decision==3);
    for j = 1:numel(si)
        s = si(j); k = ki(j);
        v = o.Cboot{s,k}; v = v(~isnan(v));
        if numel(v) < 2, continue; end

        % ecart-type parametrique sur les donnees completes, meme
        % standardisation et meme methode d'ajustement que le pipeline
        Xsk = X(:,:,s,k); [Nn,Rr] = size(Xsk);
        x = (Xsk(:)-mean(Xsk(:)))/std(Xsk(:)); y = Y(:);
        T = table(y, x, nominal(repelem((1:Nn)',Rr,1)), ...
                  'VariableNames',{'Y','X','Sample'});
        lme = fitlme(T,'Y ~ X + (1|Sample)','FitMethod','REML');

        % Ecart-type bootstrap CORRIGE par n/(n-1), tel que le pipeline le
        % rapporte (o.b1_se). Le biais normalise est calcule sur les pentes
        % brutes : la correction multipliant numerateur et denominateur par
        % le meme facteur, il n'en depend pas.
        sb(s,k) = o.b1_se(s,k);
        sp(s,k) = lme.Coefficients.SE(2);
        rt(s,k) = sb(s,k)/sp(s,k);
        bn(s,k) = (mean(v) - o.beta1(s,k))/std(v);
        sk(s,k) = skewness(v);
        kp(s,k) = true;
    end
    ratio(r,:,:)=rt; biasN(r,:,:)=bn; skw(r,:,:)=sk;
    seB(r,:,:)=sb;  seP(r,:,:)=sp;  kept(r,:,:)=kp;
end
fprintf('Termine en %.1f min.\n\n', toc(tt)/60);

% --- tableau par paire --------------------------------------------------------
freq = squeeze(mean(kept,1));
fprintf('Paires retenues dans au moins %.0f %% des replicats.\n', 100*FREQ_MIN);
fprintf('rapport = ecart-type bootstrap / ecart-type parametrique\n');
fprintf('biais   = (moyenne bootstrap - pente estimee) / ecart-type bootstrap\n\n');
fprintf(' desc   s   freq    rapport median   [q25   q75]    biais median   asym\n');
for k = 1:K
    for s = 1:S
        if freq(s,k) < FREQ_MIN, continue; end
        rv = ratio(:,s,k); rv = rv(~isnan(rv));
        bv = biasN(:,s,k); bv = bv(~isnan(bv));
        sv = skw(:,s,k);   sv = sv(~isnan(sv));
        q  = quantile(rv,[0.25 0.75]);
        fprintf(' %-5s %3d  %5.2f   %10.2f    [%4.2f  %4.2f]   %10.2f   %6.2f\n', ...
            names{k}, s, freq(s,k), median(rv), q(1), q(2), median(bv), median(sv));
    end
end

% --- stabilite du biais -------------------------------------------------------
% Un biais reproductible se corrige, un biais erratique ne se corrige pas. Le
% critere est le rapport entre la mediane du biais normalise et sa dispersion
% sur replicats : au-dessus de 1, le biais est systematique.
fprintf('\n--- Stabilite du biais sur replicats ---\n');
fprintf(' desc   s   mediane   ecart-type   |mediane|/ecart-type\n');
for k = 1:K
    for s = 1:S
        if freq(s,k) < FREQ_MIN, continue; end
        bv = biasN(:,s,k); bv = bv(~isnan(bv));
        fprintf(' %-5s %3d  %7.2f   %10.2f   %18.2f\n', ...
            names{k}, s, median(bv), std(bv), abs(median(bv))/std(bv));
    end
end

% --- figure -------------------------------------------------------------------
% Deux panneaux : rapport contre echelle, biais normalise contre echelle. Une
% courbe par descripteur, mediane et intervalle interquartile sur replicats.
kShow = [];
for k = 1:K, if any(freq(:,k) >= FREQ_MIN), kShow(end+1)=k; end, end %#ok<AGROW>
col = lines(numel(kShow));

f = figure('Position',[100 100 900 380]);
for panel = 1:2
    subplot(1,2,panel); hold on;
    for j = 1:numel(kShow)
        k = kShow(j);
        ss = find(freq(:,k) >= FREQ_MIN)';
        med = nan(size(ss)); lo = med; hi = med;
        for a = 1:numel(ss)
            if panel==1, v = ratio(:,ss(a),k); else, v = biasN(:,ss(a),k); end
            v = v(~isnan(v));
            med(a) = median(v); q = quantile(v,[0.25 0.75]);
            lo(a) = q(1); hi(a) = q(2);
        end
        fill([ss fliplr(ss)], [lo fliplr(hi)], col(j,:), ...
             'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
        plot(ss, med, '-o', 'Color',col(j,:), 'LineWidth',1.4, 'MarkerSize',5);
    end
    xlabel('Scale index s');
    if panel==1
        ylabel('se_{bootstrap} / se_{parametric}');
        yline(1,':','HandleVisibility','off');
        legend(names(kShow),'Location','best'); legend boxoff;
    else
        ylabel('Normalised bootstrap bias');
        yline(0,':','HandleVisibility','off');
        yline(-0.25,'--','HandleVisibility','off');
    end
    grid on; box on; set(gca,'FontSize',10);
end
print(f,'Fig_boot_vs_param.png','-dpng','-r300');
fprintf('\nFig_boot_vs_param.png sauvee.\n');

save('ws_run_boot_vs_param.mat','ratio','biasN','skw','seB','seP','kept', ...
     'freq','NREP','B_BOOT','FREQ_MIN','P','names');
