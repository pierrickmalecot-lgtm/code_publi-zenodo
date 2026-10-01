% version v1, 28/09/2026 : retrace de Fig_randslope_R sans relancer le calcul
%PLOT_RANDSLOPE_R  Relit ws_run_randslope_R.mat et retrace la figure avec une
%   seule etiquette d'axe x commune aux deux panneaux (les deux etiquettes
%   identiques se chevauchaient au centre).

load('ws_run_randslope_R.mat','res','GRID');

f = figure('Units','centimeters','Position',[2 2 15 7],'Color','w');
t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(GRID, res(1).pA, '-o', 'LineWidth',1.4, 'MarkerSize',5);
plot(GRID, res(2).pA, '--s','LineWidth',1.4, 'MarkerSize',5);
yline(0.5,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',GRID,'XTickLabel',string(GRID),'FontSize',9);
ylabel('P(AIC selects random slope)','FontSize',9);
legend({'sweep in R','sweep in n'},'Location','northwest','Box','off','FontSize',8);
grid on; box on; ylim([0 1]);

nexttile; hold on;
plot(GRID, res(1).dA, '-o', 'LineWidth',1.4, 'MarkerSize',5);
plot(GRID, res(2).dA, '--s','LineWidth',1.4, 'MarkerSize',5);
yline(0,':','HandleVisibility','off');
set(gca,'XScale','log','XTick',GRID,'XTickLabel',string(GRID),'FontSize',9);
ylabel('\DeltaAIC (reduced - full)','FontSize',9);
grid on; box on;

xlabel(t,'Regions per specimen R, or specimens per condition n','FontSize',9);

print(f,'Fig_randslope_R.png','-dpng','-r300');
fprintf('Fig_randslope_R.png sauvee.\n');
