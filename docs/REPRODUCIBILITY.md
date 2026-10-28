# Reproductibilité et niveaux de contrôle

## Ce qui est reproductible depuis le dépôt

Le fichier `data/compact/campaign_v04/campaign_cells_v04.csv` et le script `scripts/verify_publication_data.py` recomposent les 1 152 contrastes, les 64 scores primaires, les statistiques de Student et les résumés secondaires. Cette opération ne relance aucun tirage, épluchage, oracle ou solveur d'EDP.

Le CSV contient 8 064 lignes : `64 traductions × 2 rayons × 3 orientations × 3 contractions × 7 beta`. L'unité statistique du primaire est la traduction complète.

## Ce qui exige la grande archive

`_prepare_replication_v0_4_20261027.zip` contient les 64 fichiers `row_00_main.npz` à `row_63_main.npz` et permet les vérifications depuis les polygones sauvegardés. Cette archive de 210 150 972 octets est fournie comme supplément GitHub/Zenodo, mais n'est pas versionnée dans Git.

## Contrôle géométrique déjà réalisé

Le contrôle v0.6 détaillé porte sur huit lignes et 30 863 couches. Il ne doit pas être décrit comme un contrôle intégral des 64 trajectoires. L'implémentation QuickHull est partiellement indépendante : elle partage certaines primitives et les lignes initiales qualifiées avec le moteur principal.

## Dépendances retrouvées

Les quatre dépendances annoncées ont été retrouvées avec les SHA-256 exacts et sont copiées sous `src/recovered_dependencies/`. Les scripts historiques conservent leurs noms et leur provenance. Le mode d'emploi portable porte d'abord sur la recomposition statistique ; il ne prétend pas transformer tous les scripts historiques à chemins absolus en une nouvelle campagne autonome.

## Commande légère autorisée

```powershell
python scripts/verify_publication_data.py `
  --campaign-csv data/compact/campaign_v04/campaign_cells_v04.csv `
  --analysis-json data/compact/analysis_v04/analysis_results_v04.json `
  --output verification_recomposee.json
```

Python 3.10+, NumPy et SciPy sont requis.

## Interdiction de nouvelle campagne

La présente publication n'autorise ni nouveau tirage, ni nouvel épluchage, ni recalcul d'oracle ou d'EDP, ni extension de rayons ou de formes. Une nouvelle question scientifique nécessiterait un nouveau protocole.
