% version v8, 04/09/2026 : graines centralisees via paper_seed
% version v4, 01/08/2026 : taux d erreur mesure sur D4, panel a cinq descripteurs
%RUN_ROBUSTNESS_COMPARE  Compare les deux modes de diagnostic (residual vs pooled)
%   sur les MEMES tirages, pour documenter ce que l'option change :
%     - couverture des branches robustes (Welch / KW) par scenario
%     - taux de faux positifs (verrou 1)
%     - NOMBRE DE PAIRES RETENUES identique entre modes (message cle : le choix
%       du mode ne change pas la decision fonctionnelle, seulement le routage)
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P=paper_params(); NREP=200;
modes = {'residual','pooled'};
scen  = {'nominal','nonnormal','heteroscedastic'};

for md=1:numel(modes)
    mode=modes{md};
    fprintf('\n========== diagMode = %s ==========\n', mode);
    for c=1:numel(scen)
        sc=scen{c}; cntW=0; cntKW=0; fp=zeros(NREP,1); nRet=zeros(NREP,1);
        parfor r=1:NREP
            rs=RandStream('Threefry','Seed', paper_seed(sc, r));
            [X,Y]=gen_volet1(P,sc,rs);
            out=paper_pipeline(X,Y,P,struct('diagMode',mode));
            b=out.branch;
            % Null propre = D4 SEUL, bruit ROI. D5 porte un effet specimen :
            % il franchit le verrou 1 par pseudo-replication et est rattrape au
            % verrou 2, ce n'est donc pas un null de reference.
            g1=(out.p_fdr<P.q)&(out.eta2>=P.eta2_thr);
            fp(r)=mean(g1(:,4));
            cntW=cntW+nnz(b==1); cntKW=cntKW+nnz(b==2);
            nRet(r)=nnz(out.decision==3);
        end
        tot=NREP*P.S*P.K;
        fprintf(['%-16s : robuste %.1f%% (Welch %.1f%% + KW %.1f%%) | ' ...
            'FP=%.4f | paires retenues/rep = %.2f\n'], ...
            sc,100*(cntW+cntKW)/tot,100*cntW/tot,100*cntKW/tot,mean(fp),mean(nRet));
    end
end
fprintf(['\nLecture : comparer ligne a ligne les deux blocs. Le routage (Welch/KW)\n' ...
    'differe, mais "paires retenues/rep" doit rester quasi identique entre modes\n' ...
    '-> le choix du diagnostic ne change pas la selection fonctionnelle finale.\n']);
