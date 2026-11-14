# Orbit — FMM Simulation

Simulation galactique à deux niveaux : distribution collective et étoiles déterministes inférées, puis système orbital local. Calcul et rendu WebGPU, références CPU et tests de déterminisme.

## Démarrer

Avec Node.js 22 ou ultérieur :

```sh
npm start
```

Ouvrir <http://localhost:8080>. Le serveur écoute uniquement sur la machine locale. Le GPU demande un navigateur avec WebGPU ; les options CPU restent disponibles. Aucune compilation et aucune dépendance de production ne sont nécessaires.

## Développer

```sh
npm ci
npm test
npm run format:check
```

Le formatage concerne les modules déjà migrés et le contrôleur principal. Les références numériques et les noyaux non migrés restent inchangés. Voir [les commandes de validation](docs/VALIDATION.md), notamment les tests GPU et les mesures de performances.

## Organisation

- [Architecture actuelle et prochaines étapes](docs/ARCHITECTURE.md).
- `src/observation/` : contexte caméra, publication des images, politique de budget.
- `comparison-view.js` : orchestration et interface, en cours de découpage.
- Modules scientifiques à la racine : continu, matière persistante, formation, FMM et optique.
- `test-*.cjs`, `test-*.html`, `test-fixtures/` : tests numériques, navigateur et références.
- `scripts/` : serveur local et commande de tests.

La branche [`legacy`](https://github.com/lsarrazi/fmm-simulation/tree/legacy) conserve l'état fonctionnel antérieur à la refonte. `CONTINUUM.md`, `PLANETARY.md` et les anciens comptes rendus conservent le journal des expérimentations ; leurs passages historiques ne constituent pas tous la spécification actuelle.

Le modèle reste une approximation collective avec raffinement local, pas une intégration N-corps exacte de toute la population théorique. Les limites physiques sont décrites dans l'architecture.
