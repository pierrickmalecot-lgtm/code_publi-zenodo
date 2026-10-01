% version v8, 04/09/2026 : graines centralisees via paper_seed
%RUN_SCENARIOS  Section 4.8 : null + twoscale, 200 replicats, conventions papier.
clear; clc;
pool=gcp('nocreate'); if isempty(pool), parpool('local'); end
P=paper_params(); NREP=200;

anyPass=false(NREP,1); maxR2=zeros(NREP,1);
parfor r=1:NREP
    rs=RandStream('Threefry','Seed', paper_seed('null', r));
    [X,Y]=gen_volet1(P,'null',rs);
    out=paper_pipeline(X,Y,P);
    anyPass(r)=any(out.decision(:)==3);
    v=out.R2m(~isnan(out.R2m)); maxR2(r)=max([v;0]);
end
fprintf('NULL : prop. replicats avec >=1 paire retenue = %.4f | max R2 = %.3f\n',...
    mean(anyPass), max(maxR2));

ret=zeros(P.S,P.K); r2A=nan(NREP,1); r2B=nan(NREP,1);
parfor r=1:NREP
    rs=RandStream('Threefry','Seed', paper_seed('twoscale', r));
    [X,Y]=gen_volet1(P,'twoscale',rs);
    out=paper_pipeline(X,Y,P);
    ret=ret+(out.decision==3);
    r2A(r)=out.R2m(3,1); r2B(r)=out.R2m(8,2);
end
fprintf('TWOSCALE : R2m moyen pic A (s=3) = %.3f | pic B (s=8) = %.3f\n',...
    mean(r2A,'omitnan'), mean(r2B,'omitnan'));
disp('Freq. retention (sur 200):'); disp(round(ret/NREP,2));

% figure bi-panneau (replicats types)
rs=RandStream('Threefry','Seed', paper_seed('null',1)); [Xb,Yb]=gen_volet1(P,'null',rs);
outB=paper_pipeline(Xb,Yb,P);
rs=RandStream('Threefry','Seed', paper_seed('twoscale',1)); [Xc,Yc]=gen_volet1(P,'twoscale',rs);
outC=paper_pipeline(Xc,Yc,P);
cmap=[0.90 0.90 0.90;0.70 0.85 1.00;0.30 0.65 0.90;0.15 0.70 0.30];
dn={'D1*','D2*','D3','D4','D5','D6'};
fig=figure('Position',[80 80 1200 420],'Color','w');
tt={'Null-correlation configuration','Two-scale configuration'};
dd={outB.decision,outC.decision};
for pnl=1:2
    subplot(1,2,pnl); imagesc(dd{pnl}'); colormap(gca,cmap); caxis([-0.25 3.25]);
    colorbar('Ticks',[0.375 1.125 1.875 2.625],'TickLabels',...
        {'Excluded','Discrim. only','Discrim.+reliable','Retained'});
    set(gca,'YTick',1:6,'YTickLabel',dn,'XTick',1:P.S,'FontSize',9);
    xlabel('Scale level'); ylabel('Descriptor'); title(tt{pnl},'FontWeight','bold');
end
paper_figstyle(fig, 15, 6.5);
print(fig,'-dpng','-r300','Fig_decision_matrices_BC.png');
fprintf('Fig_decision_matrices_BC.png sauvee.\n');
