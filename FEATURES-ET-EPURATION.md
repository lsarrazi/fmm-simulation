# Référentiel central et multiplicité — mise à jour

La barre fixe contient désormais l’identité modifiable de la cible réellement suivie, son nom et les actions Suivre / Libre. Le modèle de condensation fractionne les concentrations trop étendues en domaines d’alimentation locaux ; voir PLANETARY.md pour les hypothèses et les mesures.

# Épuration appliquée — 13 novembre 2026

La proposition a été acceptée et mise en œuvre : barre fixe de simulation, panneaux Observer/Objet/Avancé, recherche commune, Nouvelle simulation, entrée principale unique. Les orbites de naissance ont été conservées à la demande de l’utilisateur. Les anciennes interfaces ont été retirées de la navigation principale ; le banc N-corps est dans validation.html, les références dans test-fixtures, les rapports dans validation-results et les sauvegardes dans archive/backups. Les moteurs scientifiques partagés restent disponibles.

La suite de ce document conserve l’inventaire et la proposition initiale pour référence ; ses formulations au futur décrivent cet état antérieur. Le modèle de naissance actualisé est documenté dans PLANETARY.md.

# Inventaire et proposition d'épuration

État inspecté le 13 novembre 2026. Ce document est une proposition : aucune suppression de fonction, d'écran ou de fichier n'a été effectuée dans cette intervention.

## Fonctions disponibles

| Ensemble | Ce qui existe | Proposition |
|---|---|---|
| Univers à deux niveaux | Galaxie continue, étoiles reconstruites localement, un système stellaire détaillé à la fois | Conserver : cœur du projet |
| Identités | Propriétés et noms déterministes ; recherche par adresse/rang ; état des rencontres pendant l'expérience | Conserver ; réunir les recherches étoile et corps en une seule commande contextuelle |
| Gravitation | Continu cinétique GPU, raffinement FMM stellaire, FMM GPU orbitale partagée, garde hors champ | Conserver ; sortir les paramètres numériques de l'écran principal |
| Échelle physique | Calibration Voie lactée, unités réelles, masses/rayons probabilistes, catégories déduites de la masse | Conserver ; documenter les approximations dans une aide |
| Population théorique | 100 milliards d'étoiles par défaut, options jusqu'à 10^24, corps locaux bornés indépendamment | Conserver ; déplacer les populations expérimentales dans Avancé, sans les supprimer |
| Systèmes | 6500–15000 corps matérialisés par défaut ; budgets fixes jusqu'à 30000 ; population inférée bien plus grande | Conserver ; afficher le budget actif seulement au niveau système |
| Caméra | Zoom continu et exact, déplacement, champ optique, hauteur, coordonnées x/y, cadrages, suivi de tous les objets | Garder zoom/déplacement/suivi/cadrage ; déplacer coordonnées exactes et hauteur dans Avancé |
| Temps | Pause/reprise, durée par seconde avec 14 unités, facteur rationnel exact, préréglages, boutons ×10/÷10, avance manuelle, format séparé d'horloge | Garder la précision et tous les ralentissements ; réunir les présentations dans un seul composant, mettre l'avance manuelle dans Avancé |
| Lumière stellaire | PSF, photosphère réelle, exposition, halo, profils d'aigrettes, angle, diffusion, largeur optique | Garder exposition et rendu ; placer les réglages de lunette dans Avancé et les masquer au niveau système où ils ne s'appliquent pas |
| Corps réfléchissants | Disques GPU, éclairage de l'hôte, face nocturne, rayon physique et contours pointillés | Amplification fictive supprimée ; conserver les contours de navigation |
| Navigation système | Repères, orbites de naissance, coronographe optionnel, retour hôte/galaxie, recherche de rang | Garder les repères et retours ; proposer de retirer les orbites de naissance de la vue normale car elles ne montrent pas les trajectoires perturbées ; coronographe dans Avancé |
| Inspecteur | Deux fiches spécialisées étoile/corps, propriétés physiques, coordonnées et vitesses, boutons de suivi et de cadrage | Un inspecteur commun avec sections adaptées ; propriétés optiques détaillées repliées |
| Calibration de visibilité | Plafond arbitraire de sources et seuil d'apparition arbitraire | Conserver ces deux commandes demandées pour la calibration ; les présenter seulement au niveau galaxie |
| Scénarios | Voie lactée, galaxie originale, nuage froid sans rotation, nuage agité sans rotation ; structure initiale | Conserver les scénarios sans rotation ; déplacer leur choix et les paramètres de départ dans Nouvelle simulation |
| Compatibilité et qualité | Moteur CPU du continu, résolutions GPU/CPU, échantillon initial, adoucissement | Réglages Avancé ; conserver pour l'instant la référence CPU et les contrôles scientifiques |
| Diagnostics | FPS, cadence réelle, nombre de corps, coût physique/rendu, masse conservée, moment cinétique, limites mémoire | Garder FPS/cadence réelle/corps ; regrouper les autres dans Diagnostics |

## Ce que je retirerais réellement après accord

1. Le doublon fmm-simulation.html : son contenu est actuellement strictement identique à index.html. Un seul point d'entrée doit ouvrir l'application actuelle, avec un accès séparé au banc de validation N-corps si on souhaite le garder.
2. L'ancienne page continuum.html et son interface propre, après vérification de ses dépendances. Elle fait doublon avec la vue continue actuelle. Ne pas supprimer continuum-solver.js : il est encore utilisé par le moteur actuel.
3. Les textes de navigation qui promettent encore une comparaison côte à côte, alors que l'interface actuelle a une seule grande vue.
4. Les commandes strictement redondantes une fois leurs remplaçantes réunies : doubles recherches d'identité, multiples présentations permanentes de la même horloge, boutons de cadrage dispersés. Garder toutes les capacités demandées, notamment le facteur exact.
5. Les anciennes images et anciens rapports de mesure dans la racine : les ranger avec les tests plutôt que les exposer comme fichiers de l'application. Conserver les résultats nécessaires à la validation.

Les sauvegardes .bak ne sont pas toutes jetables : test-shared-fmm.cjs lit encore deux anciennes versions pour comparer les calculs. Avant tout nettoyage, convertir ces références en jeux de test explicites. Une suppression aveugle casserait la validation.

## Organisation proposée de l'interface

- Barre principale : pause, cadence exacte, temps écoulé, FPS, nombre de corps.
- Panneau Observation : exposition, visibilité des étoiles au niveau galaxie ; tailles et repères au niveau système ; une recherche d'objet commune.
- Panneau Objet : nom, identité, masse, rayon, position, vitesse, suivre, cadrer, voir la surface, revenir à l'hôte.
- Nouvelle simulation : scénario, graine/configuration initiale et population théorique. Ne pas ajouter un nouveau réglage de graine avant d'en avoir besoin : il s'agit du regroupement des choix existants.
- Avancé / Diagnostics : résolution, CPU, paramètres physiques, garde, lunette et mesures détaillées.

Cette première étape épure surtout les parcours et les doublons. Elle ne retire ni la FMM partagée, ni le continu, ni les identités inférées, ni les scénarios sans rotation. Une suppression du moteur CPU ou des scénarios historiques serait une deuxième décision, avec davantage de conséquences sur les comparaisons et la compatibilité.

## Corrections réalisées séparément de l'épuration

Le nom de l'étoile hôte apparaît dans le titre du rendu, dans le repère de l'étoile et dans le bouton de retour à l'hôte.

La ceinture unique a été remplacée par une densité initiale sur une plage radiale couvrant l'intérieur, les intervalles entre les grandes masses et l'extérieur. Elle est conditionnée par leurs distances orbitales et masses : baisse partielle de probabilité près des zones de rencontre, creux à certaines résonances internes des perturbateurs massifs, faible composante coorbitale proche de ±60°, et minorité d'orbites plus excentriques. Les phases sont déterministes. Ce sont des conditions de naissance ; elles ne sont pas recalculées à chaque image et ne déplacent pas artificiellement les corps quand une planète bouge. Ensuite, la FMM fait évoluer les interactions.

Ces choix sont inspirés de structures observées : astéroïdes intérieurs et géocroiseurs, lacunes et populations coorbitales. Voir les [familles d'astéroïdes décrites par la NASA](https://science.nasa.gov/solar-system/asteroids/facts/) et les [lacunes de résonance montrées par le JPL](https://ssd.jpl.nasa.gov/diagrams/mb_hist.html). La généralisation aux systèmes inventés reste heuristique, plane et sans histoire de formation ; elle ne garantit pas la stabilité des essaims coorbitaux.

Validation : un système de test contient environ 12 % de demi-grands axes intérieurs à la première grande masse, 72 % entre la première et la dernière, 16 % au-delà. Ce sont des résultats du tirage, pas des pourcentages programmés. Le déplacement des planètes change la distribution sans changer les masses ou les noms. Les tests orbitaux, GPU et d'interface passent. Les orbites intérieures imposent des pas plus courts ; une distribution plus réaliste peut réduire la vitesse physique atteignable à budget égal.
