# Épluchage de réseau et flots de courbure à taille finie

Réplication prospective sur 64 translations appariées — Patrick Royer, chercheur indépendant non affilié.

Version publique **v0.1** — DOI : [10.5281/zenodo.23728565](https://doi.org/10.5281/zenodo.23728565).

## Résultat documenté

À `R=512`, le contraste moyen préspecifié après soustraction de la sensibilité numérique empirique favorise `beta=0.32` face à `1/3` sur 64 nouvelles translations. La borne inférieure de Student unilatérale nominale vaut `1.101658655980197e-7 > 0`. Le classement s'inverse à `R=256` et n'est pas uniforme selon l'orientation ou la contraction.

Ce résultat est une comparaison moyenne à taille finie. Il ne détermine aucun exposant exact, ne prouve aucune convergence, ne réfute pas une limite affine et ne concerne aucune nouvelle valeur de pi.

## Organisation

- `output/pdf/` : note expérimentale finale ;
- `manuscript/` : source LaTeX, bibliographie et commande de compilation ;
- `figures/` : figures utilisées dans le manuscrit ;
- `data/compact/` : données compactes sauvegardées (8 064 cellules) ;
- `scripts/` : vérificateur portable, génération de la figure primaire et contrôle du manifeste ;
- `src/frozen/` : copies des plans, manifestes et programmes gelés ;
- `src/recovered_dependencies/` : quatre dépendances retrouvées avec les empreintes attendues ;
- `evidence/v0.6/` : décision, fiche de synthèse historique v0.6 et bilans de contrôle ;
- `docs/` : portée, provenance, compilation, archives et reproductibilité ;
- `CONTROLE_PDF_FINAL.json` et `RAPPORT_CONTROLE_PDF_FINAL.md` : contrôle visuel final ;
- `MANIFEST_SHA256.txt` : empreintes du paquet versionné.

`ZENODO_METADATA_v0_1.json` documente les métadonnées du dépôt manuel. Il ne s'agit volontairement pas d'un fichier `.zenodo.json` : la publication Zenodo est effectuée une seule fois, depuis le brouillon ayant réservé le DOI, sans intégration GitHub automatique.

## Vérification statistique légère

```powershell
python scripts/verify_publication_data.py `
  --campaign-csv data/compact/campaign_v04/campaign_cells_v04.csv `
  --analysis-json data/compact/analysis_v04/analysis_results_v04.json `
  --output verification_recomposee.json
```

Ce script ne lance ni tirage, ni épluchage, ni oracle, ni EDP. Le fichier de sortie est un reçu local ; il n'est pas nécessaire de le committer.

Pour vérifier le paquet :

```powershell
python scripts/check_manifest.py
```

## Grande archive hors Git

`_prepare_replication_v0_4_20261027.zip` (210 150 972 octets, SHA-256 `f225c317be686798471e52a2f34e2bbc2fdf22406654702aeebf197797764dc3`) contient les 64 sauvegardes géométriques. Elle n'est pas dans l'historique Git ; elle est jointe à la release v0.1 et au dépôt Zenodo.

## Licences

La note, les figures et les données nouvelles de l'auteur sont diffusées sous **CC BY 4.0**. Le code original de l'auteur est diffusé sous **MIT**. Les fichiers gelés gardent leur identité et leur provenance ; aucune licence antérieure explicite n'a été retrouvée dans les paquets P06/P07/P08. La portée détaillée figure dans [LICENSE](LICENSE).
