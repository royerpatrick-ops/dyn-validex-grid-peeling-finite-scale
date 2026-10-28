# Compilation et contrôle PDF

Le manuscrit est composé avec pdfLaTeX, des polices Latin Modern incorporées, `amsmath` pour les formules et BibTeX pour les références.

```powershell
cd manuscript
./build.ps1
```

Après toute modification, y compris d'un DOI ou d'un lien :

1. nettoyer et recompiler ;
2. vérifier le journal LaTeX ;
3. rouvrir le PDF, vérifier les pages et les polices ;
4. rendre toutes les pages à 200 dpi au minimum ;
5. inspecter chaque page, puis les mathématiques et tableaux à 300 dpi si nécessaire ;
6. recalculer le SHA-256 et mettre à jour le bilan.

Le fichier `CONTROLE_PDF_FINAL.json` décrit la version effectivement examinée.

## Dépendances Python

Le vérificateur statistique portable requiert NumPy et SciPy. Matplotlib est optionnel : il n'est nécessaire que pour régénérer la figure produite par `scripts/plot_primary_scores.py`. Aucun de ces scripts ne relance la campagne, l'épluchage, l'oracle ou l'EDP.
