# Épluchage de réseau et flots de courbure à taille finie

Compendium reproductible de Patrick Royer, chercheur indépendant non affilié (Dynagénèse).

## Version 0.2 — complément π/4

La version **v0.2** du 2 novembre 2026 ajoute la note technique / prépublication complémentaire :

> **Anisotropie de réseau à taille finie dans l’épluchage convexe : résidus angulaires, certificats de sauts et prévision conditionnelle à un pas**

- note et mode d’emploi : [`supplements/pi4/`](supplements/pi4/) ;
- PDF contrôlé : [`note_pi4.pdf`](supplements/pi4/manuscript/note_pi4.pdf) ;
- release GitHub : [v0.2](https://github.com/royerpatrick-ops/dyn-validex-grid-peeling-finite-scale/releases/tag/v0.2) ;
- DOI de cette version : [10.5281/zenodo.23886806](https://doi.org/10.5281/zenodo.23886806) ;
- DOI concept : [10.5281/zenodo.23728564](https://doi.org/10.5281/zenodo.23728564).

Ce complément documente des contrôles de symétrie, la localisation signée et les compensations du mode 3, un certificat suffisant de saut unitaire et une prévision conditionnelle exacte à un pas à partir du polygone complet courant. Sur 263 364 naissances informatives, le certificat renforcé en couvre 262 661 : 258 077 complètes et 4 584 censurées ; les 703 restantes sont des abstentions. Aucun échec certifié complet n’est observé. Les 1 152 successeurs testés sont exacts ; l’erreur maximale sur `a3` dans le replay Windows est `2.1589855502161948e-16`.

Ces réanalyses sont rétrospectives et portent sur un pas microscopique. Elles n’établissent ni loi macroscopique fermée, ni temps d’atteinte, ni nouvel exposant universel, ni nouvelle valeur de π, ni preuve asymptotique générale, ni validation physique de la Dynagénèse. Les événements successifs ne sont pas supposés indépendants.

## Publication initiale v0.1

La version historique **v0.1** reste figée au commit `0d14755e1c4234a3f8004fdef4e0d89218fc7733` et au DOI [10.5281/zenodo.23728565](https://doi.org/10.5281/zenodo.23728565).

À `R=512`, le contraste moyen préspecifié après soustraction de la sensibilité numérique empirique favorise `beta=0.32` face à `1/3` sur 64 nouvelles translations. La borne inférieure de Student unilatérale nominale vaut `1.101658655980197e-7 > 0`. Le classement s’inverse à `R=256` et n’est pas uniforme selon l’orientation ou la contraction.

Ce résultat est une comparaison moyenne à taille finie. Il ne détermine aucun exposant exact, ne prouve aucune convergence, ne réfute pas une limite affine et ne concerne aucune nouvelle valeur de π.

## Organisation

- `supplements/pi4/` : complément v0.2, source LaTeX, figure, scripts portables, données compactes et contrôles ;
- `output/pdf/` et `manuscript/` : note expérimentale v0.1 et sa source ;
- `figures/` : figures de la publication initiale ;
- `data/compact/` : données compactes sauvegardées de v0.1 ;
- `scripts/` : vérificateurs portables de v0.1 ;
- `src/frozen/` et `src/recovered_dependencies/` : programmes gelés et dépendances retrouvées ;
- `evidence/` : décisions, provenance et bilans de contrôle ;
- `docs/` : portée, compilation, archives et reproductibilité ;
- `ZENODO_METADATA_v0_2.json` : instantané des métadonnées et identifiants réels de v0.2 ;
- `MANIFEST_SHA256.txt` : manifeste SHA-256 du paquet versionné v0.2.

Les cinq archives π/4 originales, qui totalisent 155 288 563 octets, ne sont pas ajoutées à l’historique Git. Elles sont distribuées comme actifs de la release v0.2 et de Zenodo. Leurs empreintes figurent dans [`supplements/pi4/README.md`](supplements/pi4/README.md).

## Vérification statistique légère de v0.1

```powershell
python scripts/verify_publication_data.py `
  --campaign-csv data/compact/campaign_v04/campaign_cells_v04.csv `
  --analysis-json data/compact/analysis_v04/analysis_results_v04.json `
  --output verification_recomposee.json
```

Ce script ne lance ni tirage, ni épluchage, ni oracle, ni EDP. Le fichier de sortie est un reçu local ; il n’est pas nécessaire de le committer.

Pour vérifier le paquet versionné :

```powershell
python scripts/check_manifest.py
```

## Licences

Les notes, textes documentaires, figures et données nouvelles de Patrick Royer sont diffusés sous **CC BY 4.0**. Son code original est diffusé sous **MIT**. Les fichiers historiques et les contenus de tiers conservent leurs propres avis ; voir [LICENSE](LICENSE).
