# Vérifications

## Référence legacy

La branche `legacy` conserve le prototype fonctionnel avant refonte, sans sauvegardes ni résultats générés. Les 36 suites Node existantes passent sur cette référence.

Les fichiers de `test-fixtures/` sont des références numériques indépendantes. Ne pas les reformater ou les mettre à jour pour faire passer un test : ils servent notamment aux comparaisons exactes de FMM et d'identités.

## Exécution

```sh
npm ci
npm test
npm run format:check
```

`npm test` découvre les `test-*.cjs` à la racine, hors lanceur GPU, et les exécute dans des processus isolés. Une erreur, un arrêt anormal ou un dépassement de durée font échouer la commande. Le workflow GitHub exécute cette suite et le contrôle de formatage.

Les tests GPU sont distincts : ils nécessitent Chromium et un environnement WebGPU fonctionnel. Le lanceur existant utilise Chrome sur Windows ; ils ne sont pas annoncés comme exécutés dans la CI Linux.

```sh
npm run test:gpu
node test-gpu-runner.cjs --publication --hardware
node test-gpu-runner.cjs --camera-transition --hardware
node test-gpu-runner.cjs --viewport --hardware
node test-gpu-runner.cjs --optics --hardware
node test-gpu-runner.cjs --stellar-system --hardware
```

## Performances

```sh
node test-gpu-runner.cjs --equiv --hardware
```

Exécuter les mesures séquentiellement sur le même GPU et le même champ, sans suite CPU concurrente. Comparer le nombre de sources, les coûts de recherche, de rendu et de présentation. Les FPS varient entre essais ; une courte mesure n'établit pas à elle seule un gain ou une absence universelle de régression.

Les profils, captures et rapports sont écrits dans `validation-results/`, ignoré par Git. Ils restent consultables localement.

## Premier lot validé — 14 novembre 2026

- 38 suites Node réussies (36 suites historiques et deux nouvelles suites pour la publication et la politique de budget).
- Contrôle de formatage réussi.
- Tests GPU réussis : navigation sans bords découverts ni calcul supplémentaire, suivi pendant calcul retardé, transitions de référentiel, dimensions natives/plein écran, correspondance optique CPU/GPU et continuité des étoiles voisines.
- Test du serveur local et des chemins des nouveaux modules réussi.
- Sur deux paires d'essais alternées du même champ : legacy 59,9 puis 59,7 FPS ; refonte 60,2 puis 58,5 FPS, environ 3655 sources. Un premier essai de refonte était plus lent (45,4 FPS) ; cette baisse ne s'est pas reproduite dans les paires alternées. Aucun gain de performance n'est revendiqué.
- Les sources des solveurs physiques, du modèle de formation, des identités persistantes et des moteurs optiques sont inchangées par rapport à legacy.
