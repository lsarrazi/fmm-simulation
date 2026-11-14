# Proposition d'optimisation du rendu stellaire

13 novembre 2026 — propositions uniquement. Le cercle rouge et l'ombre des noms sont les seules modifications visuelles appliquées dans cette intervention. Aucune optimisation du pipeline n'a été implémentée.

## Constats mesurés

Le rendu stellaire emploie déjà une seule commande instanciée pour ses sources, des halos calculés sur GPU, un noyau d'instrument mis en cache et une présentation GPU directe. Les déplacer sur GPU ou ajouter simplement de l'instanciation ne serait donc pas une nouveauté.

Le test actuel test-star-performance.html, avec 145 à 220 sources affichées dans les vues étudiées et la cadence initiale de l'application, mesure environ 8–12,5 ms par appel à PersistentMatter.sample, 0,7–1,3 ms de préparation/soumission CPU par appel au rendu et 0,06–0,1 ms pour sa présentation. L'affichage se situe autour de 54–60 FPS dans ce protocole. Les appels de physique locale se produisent moins souvent à cette cadence ; ces chiffres ne représentent pas une charge physique maximale.

Les temps du wrapper render ne sont pas le temps d'exécution complet du GPU. Le test optique séparé mesure environ 3,9 ms pour 512 étoiles sur une image de 256 pixels de côté, en attendant la fin de la file GPU. Ce test diffère de la vue principale : les valeurs ne doivent pas être additionnées ou comparées comme si elles décrivaient la même image. Rapports dans validation-results/star-performance-validation.txt et gpu-validation-optics.txt.

## Options proposées

| Priorité | Changement | Gain attendu | Précautions |
|---|---|---|---|
| 1 | Index spatial des groupes de matière : interroger les groupes susceptibles d'entrer dans la caméra au lieu de tous les parcourir à chaque observation | Réduire le coût CPU actuellement dominant de la reconstruction | Bornes conservatrices tenant compte du mouvement et des déformations ; conserver l'ordre des rangs pour ne pas modifier les identités |
| 2 | Réutiliser les tableaux de préparation/transfert et mettre en cache les propriétés immuables des identités observées, avec un budget borné | Réduire allocations, pauses de mémoire et calculs répétés ; bénéfice croissant avec le plafond de sources | Ce cache est indépendant de la sélection des sources : les règles de disparition au dézoom restent identiques |
| 3 | Borner l'emprise de chaque halo selon son flux et le seuil optique déjà utilisé | Réduire les fragments GPU parcourus pour des halos faibles dont la majorité finit invisible | Garder la normalisation du flux et vérifier les positions fractionnaires, aigrettes et transitions de zoom |
| 4, optionnelle | Séparer les noyaux/photosphères et aigrettes fins du halo large ; calculer seulement le halo diffus à résolution réduite | Diminuer le coût des grandes surfaces de halo et du mélange HDR dans les champs denses | Risque de différence visuelle ; conservation d'énergie et comparaison d'images nécessaires. Ne pas réduire la résolution des surfaces ou des aigrettes fines |

Avant de choisir la troisième ou quatrième option, compléter les mesures avec le coût GPU par passe et des plafonds plus élevés. Les premières mesures montrent le coût de préparation et d'échantillonnage ; elles ne démontrent pas encore que la rasterisation est le principal goulot sur toutes les vues.

## Recommandation

Commencer par l'index spatial et la réutilisation des données, puis mesurer de nouveau. Appliquer ensuite la limitation des emprises de halo si le profil GPU le justifie. Garder la résolution réduite des halos comme option distincte, à approuver après comparaison visuelle.

Le critère de succès serait une réduction du temps par actualisation et des pires saccades à nombre de sources, exposition, optique et scénario égaux, sans changer les identités, la distribution lumineuse ou la conservation des sources au zoom constant et au zoom avant. Aucun gain multiplicatif n'est annoncé avant mesure.
