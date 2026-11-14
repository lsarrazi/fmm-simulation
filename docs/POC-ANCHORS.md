# POC : galaxie à ancres numériques

Branche : `poc/anchor-galaxy`. Page : `poc-anchors.html`. Le moteur de `comparison.html` reste disponible et n'est pas remplacé.

## Essayer

```sh
npm ci
npm start
```

Ouvrir <http://localhost:8080/poc-anchors.html>. Un navigateur WebGPU est nécessaire. La page démarre en pause avec 10 000 ancres et 100 milliards d'identités théoriques. Choisir 4 096, 10 000 ou 16 384 ancres puis recommencer pour changer la résolution ; cela crée un autre catalogue.

Pour le diagnostic : atteindre 26,5 My, puis cliquer sur « Observer la région de la coupure ». Ce bouton cadre les coordonnées physiques du diagnostic historique, pas l'ancienne identité, qui appartient à un autre modèle. Zoomer pour découvrir les étoiles. Le mode « Densité seule » permet d'examiner la distribution sans sélection lumineuse.

Clic : inspecter ; clic molette : suivre ; glisser molette : décaler le cadrage en gardant le référentiel ; glisser gauche : caméra libre. Une identité peut être saisie dans le champ Référentiel. La vitesse par défaut est 1 année simulée par seconde ; en marche, le plafond est 100 années/s lorsque des sources sont affichées. L'avance manuelle jusqu'à une époque sert au diagnostic et travaille par lots bornés.

## Modèle

- Les ancres sont des macro-particules numériques de masse égale `M_galaxie/N`. Elles représentent toute la masse stellaire du modèle, sans être des étoiles observables exceptionnellement massives.
- Leurs conditions initiales reprennent le modèle Voie lactée existant : disque, bulbe, perturbations et dispersion de vitesse. Un halo lisse reste présent.
- La gravitation utilise les passes FMM GPU existantes. L'adoucissement est 0,03 unité modèle, soit 300 pc. L'intégration kick-drift-kick travaille avec un pas maximal de 0,1 My. Le budget d'un lot est borné à huit pas.
- Chaque ancre porte un noyau probabiliste gaussien en quatre dimensions : deux coordonnées initiales de position et deux de vitesse. Sa largeur spatiale initiale dépend de la densité, entre 80 et 1 200 pc ; sa dispersion interne initiale est 5 km/s. La dispersion attribuée aux centres est réduite pour réserver cette contribution.
- Les matrices de position `A` et de vitesse `B` évoluent avec l'équation variationnelle linéarisée autour de l'ancre. Le tenseur de marée est estimé sur 32 × 32 points cibles sans masse, via le même calcul FMM, puis différencié et interpolé.
- La densité projetée est la somme analytique des gaussiennes `A Aᵀ`. Elle n'utilise pas un flou ajouté à une distribution trouée. Les noyaux se recouvrent et n'ont pas les frontières rectangulaires strictes de l'ancien modèle.

Il n'y a aucun solveur cinétique ni tampon `(x,y,vx,vy)` discret. La petite grille de cibles échantillonne uniquement les dérivées du champ de force ; elle ne transporte pas de matière. Le calcul de densité pour l'affichage reste une rasterisation CPU bornée ; la lumière et le rendu final réutilisent le moteur optique GPU.

## Identités et recherche locale

Adresse : `anchors:v1:seed:nombreAncres:population/ancreEnBase36/rang`. Les rangs sont des BigInt. La population est répartie exactement entre les ancres, sans allouer un tableau de toutes les identités.

Un arbre binaire déterministe dans les quatre coordonnées uniformes initiales permet de retrouver un rang. Une transformation normale donne les quatre coordonnées latentes. Ces coordonnées restent inchangées : position et vitesse sont évaluées avec les matrices courantes de la même ancre. Il n'y a aucune réaffectation à l'ancre la plus proche. L'appartenance à un noyau est numérique, pas une liaison gravitationnelle à une étoile centrale.

Les quantiles extrêmes sont bornés numériquement à environ huit écarts types. La projection optique coupe les queues au-delà de six écarts types. Ces seuils ne correspondent pas à une frontière importante de densité.

La recherche intersecte des sous-arbres de rangs avec la fenêtre. Le budget de travail est réparti entre les noyaux selon leur probabilité de présence. Les cellules partiellement résolues utilisent une approximation d'intersection ; les propositions hors champ sont rejetées. L'interface signale cette approximation et le nombre de sources manquantes par rapport à l'allocation souhaitée. Le plafond lumineux n'est donc pas une garantie de complétude.

Les coordonnées de naissance ont un cache FIFO borné. Les recherches peuvent être réutilisées pendant au plus dix années simulées et tant que déplacement estimé de la matière plus déplacement de caméra restent sous 0,25 pixel CSS. Les positions et vitesses des identités sont toujours réévaluées ; les objets hors champ sont écartés. Ce cache ne conserve pas une image déplacée ni un ensemble de points figés.

## Limites à évaluer avant une adoption

Ce POC teste une représentation, sans démontrer qu'elle remplace fidèlement la dynamique du continu.

1. La déformation de chaque gaussienne reste linéaire. Elle ne décrit pas correctement tous les replis non linéaires de l'espace des phases après une évolution longue. Une extension adaptative ou un modèle de transport plus riche sera probablement nécessaire.
2. Les ancres portent la gravité adoucie ; leurs distributions gaussiennes reconstruites n'alimentent pas à leur tour une résolution de champ exacte. La cohérence entre ces deux représentations de masse est approximative.
3. La reconstruction est lissée à une résolution numérique finie. L'absence de coupure ne suffit pas à prouver une densité correcte. Il reste à comparer profils radiaux, distributions de vitesse, énergie et évolution des structures avec des références à résolution croissante.
4. La découverte 4D est actuellement incomplète à travail borné. Elle demeure le principal coût lors des déplacements ou changements d'échelle. Un noyau de densité régulier ne garantit pas à lui seul un catalogue affiché sans variation d'échantillonnage.
5. Le calcul des matrices, la découverte et la préparation de la densité utilisent encore le CPU. Le FMM et la lumière utilisent le GPU. Ce n'est pas une promesse de 60 FPS, notamment en avance rapide.
6. Pas de FMM de rencontres entre étoiles inférées, pas de systèmes planétaires et pas de migration des anciennes identités dans ce POC. Ces fonctions restent dans l'application principale.

## Validation

```sh
npm test
npm run test:anchors
npm run test:anchors-ui
node test-gpu-runner.cjs --local --hardware
```

Les tests Node couvrent les identités jusqu'à 10²⁴, l'éviction du cache, la continuité de densité, la masse projetée, la recherche locale bornée et les unités/intégration avec une référence directe indépendante. Les tests GPU utilisent 10 000 ancres jusqu'à 26,5 My, comparent le champ sur 40 cibles au calcul direct, puis vérifient l'interface et le suivi en simulation.

Premières mesures sur cette machine : écart RMS du champ d'environ 1,6 % sur les cibles testées ; dérive relative du moment cinétique d'environ 0,003 % à 26,5 My. Dans la région historique, le rapport de densité entre deux points à ±50 pc vaut environ 1,03 dans ce POC. Ce ne sont pas des comparaisons de deux univers identiques.

L'optimisation de la découverte a réduit un essai d'environ 533 à 80 ms ; une recherche réutilisée prend environ 1 ms. Le mode FMM absolu à une passe réduit l'essai d'interface d'avance à 26,5 My à environ 3,6 secondes. L'essai de suivi à 1 an/s mesurait environ 23 FPS. Ces mesures sont indicatives, distinctes des performances de l'application principale.

Les rapports complets restent dans `validation-results/`, ignoré par Git.
