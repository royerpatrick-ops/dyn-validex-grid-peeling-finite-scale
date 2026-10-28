# Compilation du manuscrit

Depuis PowerShell :

```powershell
./build.ps1
```

La commande exécute directement `pdfLaTeX`, `BibTeX`, puis deux nouveaux passages `pdfLaTeX`. Le PDF final est copié dans `../output/pdf/`.

Compilation manuelle équivalente :

```powershell
pdflatex -jobname=DYN-Note_Experimentale_Grid_Peeling_v0_1 main.tex
bibtex DYN-Note_Experimentale_Grid_Peeling_v0_1
pdflatex -jobname=DYN-Note_Experimentale_Grid_Peeling_v0_1 main.tex
pdflatex -jobname=DYN-Note_Experimentale_Grid_Peeling_v0_1 main.tex
```

Le DOI de version `10.5281/zenodo.23728565` a été réservé avant la compilation finale.
