% version v10, 02/10/2026 : objets graphiques exclus de la sauvegarde ws_*.mat
% version v9, 29/09/2026 : selection par section du papier, scripts du 27-29/09
%                          (run_bca_stability, run_r2_means_calib, run_r2_null,
%                          plot_randslope_R, make_bootstrap_table), etudes non
%                          citees isolees dans le mode archive
% version v8, 17/09/2026 : ajout de run_attenuation_c et run_bca_coverage
% version v7, 04/09/2026 : nombre de travailleurs reglable et borne au profil
% version v6, 02/09/2026 : integration des six etudes de robustesse (v7)
% version v5, 01/08/2026 : ajout du mode long
function run_all(varargin)
%RUN_ALL  Lance la chaine de calcul complete dans l'ordre.
%
%   run_all              enchaine tout, calculs longs compris
%   run_all('court')     saute les etudes Monte Carlo longues
%   run_all('long')      ne fait que ces etudes longues et leurs figures
%   run_all('figures')   ne refait que les figures, a partir des resultats
%                        deja calcules
%   run_all('sec4')      ne lance que ce qui alimente la section 4 ; de meme
%                        'sec3', 'sec5', 'sec6', 'appB', 'appC', 'appE'
%
%   La correspondance entre scripts,
%   tableaux, figures et sections est donnee dans README.md.
%
%   La console est enregistree dans un journal date, de sorte qu'une fermeture
%   de session ou un redemarrage ne detruise pas les sorties. Chaque etape est
%   chronometree et isolee : si l'une echoue, les suivantes s'executent quand
%   meme et le bilan final signale ce qui n'a pas abouti.
%
%   L'espace de travail de chaque etape est enregistre dans ws_<script>.mat des
%   qu'elle s'acheve. Les scripts commencant tous par clear, sans cette
%   sauvegarde seul le dernier laisserait une trace exploitable.
%
%   ORDRE. Le volet 1 vient en premier, il alimente le tableau par descripteur
%   et les figures principales. Viennent ensuite les scenarios complementaires,
%   qui reposent sur le meme generateur, puis les etudes de robustesse, enfin
%   les deux balayages Monte Carlo, les plus couteux.
%
%   Les six etudes ajoutees en v7 (violations d'hypotheses, courbure, contraste,
%   calibration du drapeau fragile, selection de modele) sont toutes des etudes
%   Monte Carlo a 200 replicats : elles sont classees longues et n'entrent donc
%   pas dans le mode court. make_unreliable_figure recalcule ses donnees au lieu
%   de relire un fichier, il n'est donc pas classe figure seule.

opt = 'complet';
if nargin >= 1, opt = lower(varargin{1}); end

% --- journal ------------------------------------------------------------------
stamp = datestr(now,'yyyy-mm-dd_HHMM');
logfile = sprintf('journal_%s.txt', stamp);
diary(logfile); diary on;
fprintf('=== Chaine de calcul, mode %s, %s ===\n\n', opt, datestr(now));

% --- pool parallele -----------------------------------------------------------
% NWORKERS : nombre de travailleurs souhaite. Laisser [] pour le choix
% automatique, ou fixer une valeur ici.
%
% CHOIX DE LA VALEUR. Les etudes de ce dossier enchainent des ajustements de
% modeles mixtes sur de petits tableaux : le travail est limite par le calcul,
% pas par la memoire. Deux consequences.
%   - L'hyperthreading n'apporte rien. feature('numcores') peut renvoyer le
%     nombre de THREADS et non de coeurs, 24 sur un i7-13700K par exemple, qui
%     compte 16 coeurs physiques. Demander 24 travailleurs revient a en mettre
%     deux par coeur, ce qui ralentit plutot qu'accelerer.
%   - Sur les processeurs hybrides Intel a partir de la 12e generation, les
%     coeurs E sont nettement plus lents que les coeurs P. parfor repartit les
%     iterations dynamiquement, donc les inclure reste avantageux, mais le gain
%     par travailleur supplementaire decroit.
% Point de depart raisonnable : le nombre de coeurs physiques. Sur un i7-13700K,
% 16. Si la machine sert a autre chose pendant le calcul, en retirer deux.
NWORKERS = [];

nWorkers = NWORKERS;
if isempty(nWorkers)
    nWorkers = feature('numcores');
end
try
    cl = parcluster('Processes');
    if nWorkers > cl.NumWorkers
        fprintf(['Profil Processes limite a %d travailleurs, %d demandes.\n' ...
                 'Pour lever la limite : cl = parcluster(''Processes''); ' ...
                 'cl.NumWorkers = %d; saveProfile(cl);\n'], ...
                 cl.NumWorkers, nWorkers, nWorkers);
        nWorkers = cl.NumWorkers;
    end
catch
    % profil inaccessible : on laisse parpool choisir sa valeur par defaut
    nWorkers = [];
end
pool = gcp('nocreate');
if isempty(pool)
    if isempty(nWorkers)
        fprintf('Ouverture du pool parallele, taille par defaut.\n\n');
        parpool('Processes');
    else
        fprintf('Ouverture du pool parallele, %d travailleurs.\n\n', nWorkers);
        parpool('Processes', nWorkers);
    end
elseif ~isempty(nWorkers) && pool.NumWorkers ~= nWorkers
    fprintf('Pool deja ouvert avec %d travailleurs, conserve tel quel.\n\n', ...
        pool.NumWorkers);
end

% --- definition des etapes ----------------------------------------------------
% colonne 1 : script ou commande   colonne 2 : libelle
% colonne 3 : long                 colonne 4 : figure seulement
% colonne 5 : sections du papier alimentees, separees par des espaces
%             ('4 6', 'C', ...) ; 'fig' = figures seules
% Les scripts de figures marques longs ne le sont pas en eux-memes, mais ils
% lisent les fichiers produits par les etudes Monte Carlo. Les ecarter du mode
% court evite des echecs previsibles. Une commande avec parentheses est
% evaluee telle quelle dans l'espace de travail de base.
steps = {
 'run_volet1'                , 'Jeu de reference, resultats principaux' , false, false, '4 6' ;
 'make_bootstrap_table(out)' , 'Tableau du bootstrap, paires retenues'  , false, false, '4'   ;
 'make_results_figures'      , 'Figures de resultats'                   , false, true , '4'   ;
 'run_scenarios'             , 'Scenarios null et deux echelles'        , false, false, '4'   ;
 'run_fragile'               , 'Drapeau de fragilite'                   , false, false, '4'   ;
 'run_robustness_compare'    , 'Branches sous violation d hypotheses'   , false, false, '4'   ;
 'run_value_added'           , 'Apport des verrous amont'               , false, false, '6'   ;
 'run_unreliable'            , 'Descripteur correle mais peu fiable'    , false, false, '4'   ;
 'make_unreliable_figure'    , 'Courbes du descripteur peu fiable'      , false, false, '4'   ;
 'run_loso_sweep'            , 'Calibration du drapeau fragile'         , true , false, '4'   ;
 'run_boot_vs_param'         , 'Bootstrap contre parametrique'          , true , false, '4'   ;
 'run_bca_coverage'          , 'Couverture de l intervalle de pente'    , true , false, '4'   ;
 'run_bca_stability'         , 'Stabilite de l intervalle BCa'          , true , false, '4'   ;
 'run_perm_compare'          , 'Branchement contre permutation'         , true , false, '4'   ;
 'run_perm_seeds'            , 'Stabilite du niveau a n = 5'            , true , false, '4'   ;
 'run_attenuation'           , 'Correction d attenuation, cas (b)'      , true , false, '3 5' ;
 'run_attenuation_c'         , 'Correction d attenuation, cas (c)'      , true , false, '3 5' ;
 'run_latent'                , 'z0 hors du cas exact'                   , true , false, '3'   ;
 'run_stress_tests'          , 'Violations d hypotheses'                , true , false, '5'   ;
 'run_perm_size'             , 'Taille d echantillon et permutation'    , true , false, '5'   ;
 'run_curvature_noise'       , 'Courbure croisee avec le bruit'         , true , false, '5'   ;
 'run_ushape'                , 'Reponse a optimum interieur'            , true , false, '5'   ;
 'run_volet2_pipeline'       , 'Dimensionnement experimental'           , true , false, '5 C' ;
 'make_sensitivity_figures'  , 'Figures de dimensionnement'             , true , true , '5 C' ;
 'run_contrast_sweep'        , 'Contraste et grille de fiabilite'       , true , false, '6'   ;
 'run_r2_null'               , 'R2 obtenu par hasard, sans lien X-Y'    , true , false, '6'   ;
 'run_r2_means_calib'        , 'Seuil R2 du protocole simplifie'        , true , false, 'B'   ;
 'run_randslope'             , 'Selection intercept ou pente aleatoire' , true , false, 'E'   ;
 'run_randslope_R'           , 'Puissance sur la pente aleatoire'       , true , false, 'E'   ;
 'make_all_figures'          , 'Retrace des figures depuis ws_*.mat'    , true , true , 'fig' ;
 'plot_randslope_R'          , 'Figure de la pente aleatoire'           , true , true , 'E'   ;
};

% --- selection selon le mode --------------------------------------------------
tags    = steps(:,5);
estLong = cell2mat(steps(:,3));
estFig  = cell2mat(steps(:,4));
archive = strcmp(tags,'x');
aTag    = @(t) cellfun(@(c) any(strcmp(strsplit(c,' '), t)), tags);
switch opt
    case 'complet'
        keep = ~archive;
    case 'court'
        keep = ~estLong & ~archive;
    case 'long'
        keep = estLong & ~archive;
    case 'figures'
        keep = estFig & ~archive;
    case {'sec3','sec4','sec5','sec6'}
        keep = aTag(opt(4));
    case {'appb','appc','appe'}
        keep = aTag(upper(opt(4)));
    otherwise
        error(['Mode inconnu : %s. Attendu complet, court, long, figures, ' ...
               'sec3 a sec6, appB, appC ou appE.'], opt);
end
steps = steps(keep,:);
if isempty(steps)
    error('Aucune etape pour le mode %s.', opt);
end

% --- execution ----------------------------------------------------------------
n = size(steps,1);
etat = cell(n,1); duree = zeros(n,1);
tGlobal = tic;

for i = 1:n
    nom = steps{i,1}; lib = steps{i,2};
    fprintf('\n');
    fprintf('%s\n', repmat('=',1,72));
    fprintf('  Etape %d sur %d : %s\n', i, n, lib);
    fprintf('  Script : %s\n', nom);
    fprintf('%s\n\n', repmat('=',1,72));

    estCmd = contains(nom,'(');
    base   = strtok(nom,'(');
    if exist(base,'file') ~= 2
        fprintf('  ABSENT : le fichier %s.m est introuvable.\n', base);
        etat{i} = 'absent'; continue;
    end

    t0 = tic;
    try
        % Les scripts sont lances dans l'espace de travail de base. Plusieurs
        % commencent par clear, ce qui effacerait les variables de cette
        % fonction s'ils y etaient executes directement.
        if estCmd
            evalin('base', nom);
        else
            evalin('base', sprintf('run(''%s'')', nom));
        end
        duree(i) = toc(t0);
        etat{i} = 'ok';
        fprintf('\n  Termine en %s.\n', dureeTexte(duree(i)));

        % Sauvegarde de l'espace de travail de l'etape. Chaque script commence
        % par clear, donc seul le contenu de l'etape qui vient de s'achever est
        % enregistre. Sans cela, seul le dernier script laisse une trace et les
        % resultats intermediaires sont perdus.
        % Les objets graphiques (figures, axes) sont exclus : les enregistrer
        % alourdit le fichier sans rien apporter, les figures etant deja
        % exportees en PNG par les scripts.
        wsfile = sprintf('ws_%s.mat', base);
        try
            w = evalin('base','whos');
            cl = {w.class};
            ok = ~(startsWith(cl,'matlab.ui') | startsWith(cl,'matlab.graphics'));
            if any(ok)
                lst = sprintf(',''%s''', w(ok).name);
                evalin('base', sprintf('save(''%s''%s)', wsfile, lst));
            end
            fprintf('  Espace de travail enregistre dans %s\n', wsfile);
        catch
            fprintf('  Espace de travail non enregistre.\n');
        end
    catch ME
        duree(i) = toc(t0);
        etat{i} = 'echec';
        fprintf('\n  ECHEC apres %s : %s\n', dureeTexte(duree(i)), ME.message);
        if ~isempty(ME.stack)
            fprintf('  Dans %s a la ligne %d.\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

% --- bilan --------------------------------------------------------------------
fprintf('\n\n%s\n', repmat('=',1,72));
fprintf('  Bilan, duree totale %s\n', dureeTexte(toc(tGlobal)));
fprintf('%s\n', repmat('=',1,72));
for i = 1:n
    fprintf('  %-8s %-38s %s\n', etat{i}, steps{i,2}, dureeTexte(duree(i)));
end
fprintf('\n  Espaces de travail : ws_<script>.mat, un par etape aboutie.\n');
nEchec = sum(strcmp(etat,'echec')) + sum(strcmp(etat,'absent'));
if nEchec == 0
    fprintf('\n  Toutes les etapes ont abouti.\n');
else
    fprintf('\n  %d etape(s) sans resultat, voir le detail ci-dessus.\n', nEchec);
end
fprintf('  Journal enregistre dans %s\n', logfile);
diary off;
end

% -----------------------------------------------------------------------------
function s = dureeTexte(t)
if t < 60
    s = sprintf('%.0f s', t);
elseif t < 3600
    s = sprintf('%.1f min', t/60);
else
    s = sprintf('%.1f h', t/3600);
end
end
