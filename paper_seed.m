% version v2, 17/09/2026 : bloc 440-459 pour run_attenuation_c, run_bca_coverage au niveau 0
% version v1, 04/09/2026 : attribution centralisee des graines
function s = paper_seed(scenario, r, level)
%PAPER_SEED  Graine du generateur, determinee par le scenario et non par le
%   script appelant.
%
%   POURQUOI. Jusqu'a la version v7 chaque script fixait sa propre base de
%   graines. Deux consequences, toutes deux corrigees ici.
%
%   1. Le meme scenario simule sous deux scripts differents produisait des
%      jeux differents, donc des moyennes sur replicats legerement
%      differentes : le nombre moyen de paires retenues sous 'nominal'
%      valait 8.82, 8.95 ou 8.85 selon le script consulte. Rien n'etait faux,
%      mais le papier citait trois valeurs pour une meme quantite.
%
%   2. Les bases choisies script par script se recouvraient. Mesure sur la
%      v7 : run_value_added et run_loso_sweep partageaient 100 graines,
%      run_boot_vs_param et run_contrast_sweep 100 autres, et ces deux-la
%      tournent tous deux sous 'nominal', donc simulaient les memes donnees
%      sans que rien ne l'indique ; run_threshold_sensitivity et
%      run_curvature_noise en partageaient 190.
%
%   REGLE. s = 10^6 * identifiant_scenario + 10^3 * niveau + replicat.
%   La graine ne depend donc que de ce qui est simule, jamais de l'endroit
%   d'ou l'appel est fait. Deux scripts qui simulent le meme scenario au meme
%   niveau de balayage et au meme numero de replicat obtiennent le meme jeu,
%   ce qui rend leurs resultats directement comparables. Deux scenarios, ou
%   deux niveaux de balayage, sont toujours disjoints.
%
%   ARGUMENTS.
%     scenario : nom du scenario, tel qu'il est passe a gen_volet1, ou
%                'fragile' pour le generateur d'anomalie de gen_fragile.
%     r        : numero de replicat, 1 a 999.
%     level    : indice de balayage, 0 a 999. Facultatif, 0 par defaut.
%
%   ALLOCATION DES NIVEAUX. Le niveau 0 est celui des etudes SANS balayage :
%   elles partagent donc toutes le meme jeu de donnees pour un scenario et un
%   replicat donnes, ce qui est l'effet recherche. Les etudes a balayage
%   recoivent chacune un bloc disjoint, de sorte qu'un meme indice de niveau
%   ne designe jamais deux balayages differents.
%
%     0        etudes sans balayage (run_volet1, run_scenarios,
%              run_value_added, run_q_scale_selection, run_boot_vs_param,
%              run_robustness_compare, run_perm_compare, run_unreliable,
%              make_unreliable_figure, run_ushape volets 1 et 2,
%              run_bca_coverage)
%     1-19     run_stress_tests, indice de configuration
%     20-39    run_contrast_sweep, balayage du contraste
%     40-59    run_contrast_sweep, balayage du nombre de conditions
%     60-79    run_randslope, balayage de sigma_u1
%     80-99    run_loso_sweep, balayage de l'amplitude d'anomalie
%     100-119  run_perm_size, balayage du nombre de specimens
%     120-139  run_perm_seeds, jeux de graines independants
%     140-199  run_curvature_noise, balayage croise courbure x bruit
%     200-219  run_threshold_sensitivity, niveaux de bruit
%     230-239  run_randslope_R, balayage du nombre de regions
%     240-249  run_randslope_R, balayage du nombre de specimens
%     260-269  run_attenuation, balayage de la dispersion entre regions
%     270-279  run_latent, balayage du bruit latent
%     300-399  run_volet2_pipeline, configurations de dimensionnement
%     400-419  run_ushape, balayage de la position de l'optimum
%     420-439  run_ushape, volet du reajustement quadratique
%     440-459  run_attenuation_c, plans equilibre (440+) et desequilibre (450+)
%     460-504  run_r2_means_calib, pente vraie x nombre de regions
%     510-519  run_r2_null, nombre de specimens sans lien X-Y
%     460-504  run_r2_means_calib, pente vraie x nombre de regions
%     510-519  run_r2_null, nombre de specimens sans lien X-Y
%
%   NON CONCERNE. La graine 40007 de run_volet1 et run_unreliable designe le
%   TIRAGE DE REFERENCE, celui de l'encadre et des figures principales. Elle
%   est conservee telle quelle : la changer modifierait toutes les valeurs
%   publiees de l'exemple travaille sans rien apporter.
%   Les flux de reechantillonnage bootstrap et de permutation ne sont pas
%   concernes non plus : ils ne generent pas de donnees, ils reechantillonnent
%   des donnees deja generees.
if nargin < 3 || isempty(level), level = 0; end

ids = { 'nominal',         1 ; ...
        'null',            2 ; ...
        'twoscale',        3 ; ...
        'unreliable',      4 ; ...
        'nonnormal',       5 ; ...
        'heteroscedastic', 6 ; ...
        'nonlinear',       7 ; ...
        'autocorr',        8 ; ...
        'curvature',       9 ; ...
        'randslope',      10 ; ...
        'ushape',         11 ; ...
        'fragile',        12 };

k = find(strcmpi(scenario, ids(:,1)), 1);
if isempty(k)
    error('paper_seed:scenario', ...
          'Scenario inconnu : %s. Ajouter son identifiant dans paper_seed.', scenario);
end
if r < 1 || r > 999
    error('paper_seed:replicate', 'Numero de replicat hors plage 1-999 : %d', r);
end
if level < 0 || level > 999
    error('paper_seed:level', 'Indice de balayage hors plage 0-999 : %d', level);
end

s = 1e6*ids{k,2} + 1e3*level + r;
end
