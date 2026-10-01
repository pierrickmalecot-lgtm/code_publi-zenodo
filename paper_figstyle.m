% version v1, 05/09/2026 : mise en forme centralisee des figures
function paper_figstyle(fig, largeur_cm, hauteur_cm)
%PAPER_FIGSTYLE  Dimensionne une figure a sa taille d'impression et harmonise
%   les polices.
%
%   POURQUOI. Une figure MATLAB est dessinee sur un canevas dont la taille est
%   independante de celle a laquelle elle sera imprimee. Les polices sont
%   exprimees en points sur ce canevas : si le canevas fait 39 cm et que la
%   figure est inseree sur 15 cm, une police de 9 pt devient 3.5 pt sur la page.
%   Mesure sur la version v7 du dossier, aucune figure n'atteignait la taille de
%   texte demandee par l'editeur, la plupart tombant entre 3.5 et 6.3 pt.
%
%   Augmenter les polices sans toucher au canevas ne resout rien : il faudrait
%   declarer 20 a 25 pt, et le rapport entre le texte et les traits resterait
%   faux. La seule mise en forme correcte est de dessiner la figure a sa
%   largeur d'impression.
%
%   CONSIGNE SUIVIE. IOP demande un texte de 8 a 12 pt A LA TAILLE FINALE de la
%   figure, soit environ 8.5 cm pour une figure sur une colonne et 15 cm pour
%   une figure sur deux colonnes. Toutes les figures de ce dossier sont
%   composees pour 15 cm, avec une police de base de 10 pt, soit le milieu de
%   la fourchette demandee.
%
%   ARGUMENTS.
%     fig          handle de figure.
%     largeur_cm   largeur d'impression visee. Defaut 15.
%     hauteur_cm   hauteur. Defaut : la moitie de la largeur, proportion
%                  habituelle d'un graphe a un ou deux panneaux. A preciser
%                  pour les figures a plusieurs rangees.
%
%   USAGE. Appeler juste avant print, une fois tous les axes crees :
%       paper_figstyle(f, 15, 7);
%       print(f,'Fig_xxx.png','-dpng','-r300');
%
%   L'inclusion dans le document doit alors se faire SANS mise a l'echelle,
%   c'est-a-dire width=15cm et non width=\textwidth, faute de quoi le travail
%   fait ici est annule par LaTeX.
if nargin < 2 || isempty(largeur_cm), largeur_cm = 15; end
if nargin < 3 || isempty(hauteur_cm), hauteur_cm = largeur_cm/2; end

FS_BASE  = 10;   % axes, graduations
FS_LABEL = 10;   % etiquettes d'axes
FS_TITLE = 10;   % titres de panneaux
FS_LEG   = 9;    % legendes, un point sous le reste
FS_TEXT  = 9;    % annotations posees dans les axes

set(fig, 'Units','centimeters', ...
         'Position',[2 2 largeur_cm hauteur_cm], ...
         'Color','w', ...
         'PaperUnits','centimeters', ...
         'PaperPosition',[0 0 largeur_cm hauteur_cm], ...
         'PaperSize',[largeur_cm hauteur_cm]);

ax = findall(fig,'Type','axes');
for a = 1:numel(ax)
    set(ax(a),'FontSize',FS_BASE,'LineWidth',0.6,'Box','on');
    set(get(ax(a),'XLabel'),'FontSize',FS_LABEL);
    set(get(ax(a),'YLabel'),'FontSize',FS_LABEL);
    set(get(ax(a),'ZLabel'),'FontSize',FS_LABEL);
    set(get(ax(a),'Title'), 'FontSize',FS_TITLE,'FontWeight','bold');
end

lg = findall(fig,'Type','legend');
for a = 1:numel(lg)
    set(lg(a),'FontSize',FS_LEG,'Box','off');
end

cb = findall(fig,'Type','colorbar');
for a = 1:numel(cb)
    set(cb(a),'FontSize',FS_BASE);
    set(get(cb(a),'Label'),'FontSize',FS_LABEL);
end

% Textes poses dans les axes : valeurs annotees, etiquettes de points. Ils sont
% ramenes a FS_TEXT sans distinction, une figure composee a 15 cm n'ayant pas
% la place d'en accueillir de plus gros.
tx = findall(fig,'Type','text');
for a = 1:numel(tx)
    if isempty(get(tx(a),'Tag'))
        set(tx(a),'FontSize',FS_TEXT);
    end
end
end
