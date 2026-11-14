# Architecture actuelle

Orbit est une application JavaScript de navigateur. WebGPU prend en charge le continu cinétique, les interactions gravitationnelles locales et le rendu optique. Node.js sert aux tests et au serveur de développement. Prettier est une dépendance de développement uniquement.

## Chemin d'une image

1. Le solveur collectif avance la distribution en `(x, y, vx, vy)`.
2. `persistent-matter.js` retrouve les identités et positions dans la région observée ; le raffinement local ajoute les corrections de rencontre.
3. L'observation sélectionne les sources selon la caméra, leur luminosité et le budget. La garde physique reste indépendante du budget lumineux.
4. Le moteur optique produit une image CPU ou GPU.
5. `FramePublication.publish()` remplace atomiquement la description de la dernière image complète.
6. La présentation utilise cette description pour la caméra, les repères et les interactions. Un déplacement ou un dézoom qui découvrirait une région absente attend la prochaine image. Le zoom avant peut reprojeter l'image disponible lorsqu'elle couvre le cadre.

La boucle de calcul et la boucle de présentation restent distinctes. La publication n'ajoute ni rendu, ni lecture GPU, ni suréchantillonnage spatial.

## Responsabilités

| Ensemble              | Sources                                                                           | Responsabilité                                                              |
| --------------------- | --------------------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Application           | `comparison-view.js`, `comparison.html`, `comparison-layout.css`                  | Orchestration et interface, encore partiellement couplées                   |
| Caméra et publication | `src/observation/camera-frame.js`                                                 | Contexte de dimensions/unités, transformations, couverture et image publiée |
| Politique stellaire   | `src/observation/stellar-budget.js`                                               | Calcul pur du plafond adaptatif, indépendant du DOM                         |
| Horloge et outils     | `observation-tools.js`                                                            | Temps exact, visibilité, zoom ; réexport de compatibilité des outils caméra |
| Galaxie               | `continuum-solver.js`, `continuum-gpu.js`, `galaxy-comparison.js`, `milky-way.js` | Distribution collective et conditions initiales                             |
| Objets inférés        | `persistent-matter.js`, `stellar-population.js`                                   | Adresses et propriétés déterministes, recherche bornée                      |
| Gravitation locale    | `local-stellar-dynamics.js`, `fmm-core.js`, `fmm-gpu-compute.js`                  | Raffinement stellaire et FMM partagée                                       |
| Systèmes              | `formation-model.js`, `planetary-systems.js`                                      | Naissance probabiliste et système orbital actif                             |
| Optique               | `stellar-light.js`, `stellar-light-gpu.js`                                        | Photométrie, surfaces, PSF, composition HDR                                 |

## Propriété des données publiées

Une image publiée regroupe `camera`, `context`, `points`, `projection`, `solarScene` et `revision`.

- La caméra et son point d'ancrage sont copiés puis figés. Les changements de navigation ne peuvent pas les modifier.
- Le contexte de rendu est déjà immuable.
- Les tableaux d'objets et de densité sont empruntés sans copie profonde. Le producteur doit les remplacer plutôt que les modifier tant qu'une image les utilise. Le gel du conteneur ne rend pas ces tableaux immuables.
- Une réinitialisation efface la publication. Les protections de génération existantes continuent d'invalider les calculs asynchrones obsolètes.

Le GPU conserve son propre stockage d'image. `FramePublication` décrit l'image courante ; il ne constitue pas une archive d'images GPU.

## Unités et frontières

Le continu utilise les unités du modèle galactique ; la calibration les relie aux parsecs. Les systèmes utilisent UA, années et masses solaires. L'horloge conserve un temps exact en nanosecondes avec BigInt. Le contexte distingue pixels CSS, pixels natifs et largeur physique du champ.

Les fonctions d'observation pures reçoivent des paramètres validés. Les commandes et le DOM restent dans l'application. Les solveurs ne doivent pas dépendre du panneau de réglages. La catégorie affichée d'un corps ne doit pas modifier son identité ou ses équations.

## Limites physiques conservées

- Simulation dynamique plane ; la profondeur utilisée pour l'observation ne constitue pas une intégration gravitationnelle 3D.
- Les trajectoires inférées sont une approximation collective avec raffinement local borné, pas un N-corps global exact.
- Un seul système orbital actif. Les périodes inactives utilisent une propagation approchée, sans intégrer toutes les rencontres invisibles.
- La matière orbitale non matérialisée n'ajoute pas actuellement son propre champ au solveur du système.
- Le modèle de formation et les populations sont synthétiques ; leur calibration observationnelle reste partielle.

## Migration progressive

Le premier lot extrait la publication et le calcul du budget, rend le contrôleur lisible et automatise les vérifications. Les modules anciens restent à la racine pour éviter de changer simultanément tous les chemins et les conventions d'exécution. Le navigateur conserve temporairement des scripts ordonnés ; les consommateurs Node gardent les exports CommonJS historiques.

Prochains lots, chacun validé séparément :

1. Extraire navigation et référentiel commun, puis inspecteur.
2. Construire l'observation complète à une frontière explicite ; séparer davantage l'interface de l'orchestration physique.
3. Centraliser les conversions et la configuration.
4. Migrer les imports et les dossiers restants ; classer la documentation historique et automatiser les versions des ressources.

Les réécritures des noyaux numériques, les nouvelles abstractions de moteur et les changements de physique ne font pas partie de cette migration.
