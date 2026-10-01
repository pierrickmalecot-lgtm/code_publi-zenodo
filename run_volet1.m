%RUN_VOLET1  Volet 1 conforme au papier : nominal + figures + encadre 4.2 + LOSO.
clear; close all; clc;
rng(42,'twister');                       % reproductibilite (comme volet 1)
   rs=RandStream('Threefry','Seed',40007); 
P = paper_params();
[X,Y]=gen_volet1(P,'nominal',rs);      % [] -> rng global (fidele au volet 1)
out = paper_pipeline(X,Y,P,struct('bootstrap',true,'B',1000,'loso',true));

descriptorNames={'D1*','D2*','D3','D4','D5','D6'};

% ---- Encadre 4.2 (D1*, s=6) ----
s=P.s_star; k=1; m=P.n*P.R; G=P.G;
sampleCondition=repelem(1:G,P.n);
Xp=zeros(G,m); for g=1:G, Xp(g,:)=reshape(X(sampleCondition==g,:,s,k),1,m); end
xg=mean(Xp,2); xx=mean(Xp(:));
MSb=m/(G-1)*sum((xg-xx).^2); MSw=sum(sum((Xp-xg).^2))/(G*(m-1));
fprintf('\n=== ENCADRE 4.2 (D1*, s=6) — conventions papier ===\n');
fprintf('condition means : %s | global mean = %.3f\n', mat2str(round(xg',3)), xx);
fprintf('MSb = %.3f | MSw = %.4f | F = %.1f | p_FDR %s | eta2 = %.3f\n', ...
    MSb, MSw, out.F(s,k), fmtp(out.p_fdr(s,k)), out.eta2(s,k));
if abs(out.Ftest(s,k)-out.F(s,k))>1e-9
    fprintf('  (p-value issue de la branche robuste : statistique testee = %.1f)\n', out.Ftest(s,k));
end
fprintf('ICC = (F-1)/(F+%d) = %.3f | R2 marginal = %.3f | C = %.3f\n', ...
    m-1, out.ICC(s,k), out.R2m(s,k), out.C(s,k));
fprintf('LOSO : %d paire(s) fragile(s) parmi les retenues.\n', ...
    nnz(out.fragile & out.decision==3));

% ---- Figure 6 panneaux (charte volet 1) ----
maps={log10(max(out.F,10^-0.5)),'(a) log_{10}(F-statistic)',[];
      out.ICC,'(b) ICC (reliability)',[0 1];
      max(out.R2m,0),'(c) R^2 marginal (LME)',[0 1];
      out.rho,'(d) Spearman \rho',[-1 1];
      out.C,'(e) Bootstrap consistency',[0 1];
      out.eta2,'(f) \eta^2 (effect size)',[0 1]};
fig1=figure('Units','centimeters','Position',[2 2 36 20],'Color','w');
for pnl=1:6
    subplot(2,3,pnl);
    D=maps{pnl,1}; D=D'; 
    if pnl==5, nanM=isnan(D); D(nanM)=-0.05; imagesc(D,[-0.05 1]);
        colormap(gca,[0.75 0.75 0.75; parula(256)]);
    else
        if isempty(maps{pnl,3}), imagesc(D); else, imagesc(D,maps{pnl,3}); end
        colormap(gca,parula(256));
    end
    colorbar; set(gca,'YTick',1:6,'YTickLabel',descriptorNames,'FontSize',9);
    xlabel('Scale level'); ylabel('Descriptor');
    title(maps{pnl,2},'FontWeight','bold');
end
sgtitle('Multi-scale framework — Synthetic validation (paper-conformant pipeline)', ...
    'FontSize',13,'FontWeight','bold');
print(fig1,'-dpng','-r600','Fig_synthetic_example_results_hq.png');

% ---- Matrice de decision (charte volet 1) ----
fig2=figure('Units','centimeters','Position',[2 2 20 10],'Color','w');
imagesc(out.decision');
colormap(gca,[0.90 0.90 0.90;0.70 0.85 1.00;0.30 0.65 0.90;0.15 0.70 0.30]);
caxis([-0.25 3.25]);
colorbar('Ticks',[0.375 1.125 1.875 2.625],'TickLabels', ...
 {'Excluded','Discrim. only','Discrim. + reliable','Retained (all criteria)'});
set(gca,'YTick',1:6,'YTickLabel',descriptorNames,'XTick',1:P.S,'FontSize',9);
xlabel('Scale level'); ylabel('Descriptor');
title('Three-stage filtering decision matrix','FontSize',12,'FontWeight','bold');
print(fig2,'-dpng','-r600','Fig_decision_matrix.png');

% ---- Q indicator (5.8) ----
Q=zeros(P.S,6);
for s2=1:P.S, for k2=1:5
    qd=min(out.eta2(s2,k2)/0.14,1); qi=min(out.ICC(s2,k2)/P.icc_thr,1);
    if isnan(out.C(s2,k2)), qb=0; else, qb=min(out.C(s2,k2)/P.C_thr,1); end
    Q(s2,k2)=(qd+qi+qb)/3;
end, end
ret=out.decision==3;
if any(ret(:)), fprintf('Q global (paires retenues) = %.3f\n', mean(Q(ret))); end
fig3=figure('Units','centimeters','Position',[2 2 18 10],'Color','w');
imagesc(Q'); colorbar; caxis([0 1]); colormap(gca,parula);
set(gca,'YTick',1:6,'YTickLabel',descriptorNames,'XTick',1:P.S,'FontSize',9);
xlabel('Scale level'); ylabel('Descriptor');
title('Composite reliability indicator Q_{s,k}','FontSize',12,'FontWeight','bold');
print(fig3,'-dpng','-r600','Fig_Q_indicator.png');
fprintf('Figures sauvees.\n');

function s = fmtp(p)
if p==0 || p<1e-15, s='< 1e-15'; elseif p<1e-3, s=sprintf('< 1e-%d',ceil(-log10(p)));
else, s=sprintf('= %.3f',p); end
end
