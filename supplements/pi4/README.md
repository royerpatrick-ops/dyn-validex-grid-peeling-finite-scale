# Anisotropie de réseau à taille finie dans l’épluchage convexe : résidus angulaires, certificats de sauts et prévision conditionnelle à un pas

**Patrick Royer — chercheur indépendant non affilié (Dynagénèse) — 2 novembre 2026**

Note technique / prépublication complémentaire, en français avec résumé anglais.

## Liens publics

- DOI de cette version : [10.5281/zenodo.23886806](https://doi.org/10.5281/zenodo.23886806)
- DOI concept : [10.5281/zenodo.23728564](https://doi.org/10.5281/zenodo.23728564)
- publication initiale v0.1 : [10.5281/zenodo.23728565](https://doi.org/10.5281/zenodo.23728565)
- release GitHub : [v0.2](https://github.com/royerpatrick-ops/dyn-validex-grid-peeling-finite-scale/releases/tag/v0.2)
- PDF de la release : [DYN-Note_Complement_Pi4_v0_2.pdf](https://github.com/royerpatrick-ops/dyn-validex-grid-peeling-finite-scale/releases/download/v0.2/DYN-Note_Complement_Pi4_v0_2.pdf)

## Portée

Le complément étudie les effets de réseau à taille finie dans les trajectoires d’épluchage convexe déjà publiées. Il documente les contrôles de symétrie, la localisation du résidu signé du mode 3, un certificat suffisant de saut unitaire avec abstentions explicites et le transfert angulaire exact d’une prévision conditionnelle à un pas.

Sur 263 364 naissances informatives, le certificat renforcé `M > r0 + 2(L - 1)`, avec `L = ceil(m/2)`, en couvre 262 661 : 258 077 complètes et 4 584 censurées ; les 703 restantes sont des abstentions. Aucun échec certifié complet n’est observé et tous les cas suivis à π/4 sont certifiés. La preuve porte sur les points réellement présents ; la saturation permet de déduire les populations `m`, `M`, `r0` et le prochain état depuis le seul polygone courant.

Les 1 152 successeurs globaux sont reconstruits exactement à partir des polygones complets courants. L’erreur maximale `a3` du replay Windows est `2.1589855502161948e-16`. Cet accord mesure la reproduction des sorties sauvegardées, non une précision physique ni une certification de l’EDP de référence.

Les états pré-cibles sont sélectionnés rétrospectivement. Le résultat ne constitue ni une loi macroscopique autonome depuis la forme initiale, ni un temps d’atteinte. Les événements successifs ne sont pas statistiquement indépendants. Les deux reconstructions clairsemées testées à seize directions échouent à l’échelle du signal ; cela n’exclut pas toute autre fermeture.

À `π/4, R=512`, les résidus moyens sont :

| Cible | Résidu moyen |
|---|---:|
| `q=.95` | `+1.415869630128e-5` |
| `q=.90` | `-2.311500019412e-4` |
| `q=.85` | `-2.822584133007e-4` |

Les « 99 % » de localisation portent sur un petit solde net et non sur une magnitude absolue. Les fenêtres ne représentent que 5,33 % de la contribution nette du dernier pas interpolé à `q=.95` ; aucune résonance n’est établie.

## Contenu versionné dans Git

- `manuscript/note_pi4.pdf` : PDF autonome de huit pages ;
- `manuscript/note_pi4.tex` : source LaTeX avec bibliographie intégrée ;
- `manuscript/figure_prediction_pi4.png` : figure incluse dans la note ;
- `reproducibility/` : scripts portables et jeux d’entrée compacts nécessaires aux replays documentés ;
- `verification/CONTROLE_PDF_FINAL.json` : contrôle visuel du PDF exact ;
- `verification/reference/` : reçus de référence pour les certificats, la prévision et le diagnostic à seize directions ;
- `verification/windows/` : manifestes et rapprochements des résultats Windows reçus ;
- `verification/MANIFEST_SHA256_PUBLIC.json` : empreintes de la sélection publique du complément.

Les scripts historiques à chemins absolus restent dans les archives comme traces de provenance. Ils ne sont pas présentés comme exécutables autonomes.

## Archives distribuées hors de l’historique Git

| Archive | Octets | SHA-256 |
|---|---:|---|
| `DYN-Controle_Symetries_Pi4_v0_1.zip` | 18 945 079 | `028b63a97a8b01187a04be840a4ca778aa33b497114ed7fd6b230140da520242` |
| `DYN-Localisation_Angulaire_Mode3_Pi4_v0_1.zip` | 19 507 630 | `ddb3429d4fea4624b14d1a9cfacebcfcb1bab64b2c918a48dc78146ebca085f1` |
| `DYN-Deplacements_Repere_Fixe_Pi4_v0_1.zip` | 66 232 366 | `0898cb5a9cc6ceb52929bbaeec7db24ec72320acc9be0ffdfd2869635c0f566f` |
| `DYN-Sauts_Mode3_Pi4_Prediction_Conditionnelle_v0_1.zip` | 25 709 536 | `d40985bfb59fce13441a049eb2a882cf8ce82deb9be5be9be70238a8b347099e` |
| `DYN-Pi4_Resultats_Windows_Verifies_20261102.zip` | 24 893 952 | `ea1364e07c00110ad898e924bc7a82eb5a801c3983c996a5a2bfb3deedbc6956` |

Les cinq ZIP ont passé un contrôle CRC intégral. Les quatre archives scientifiques restent inchangées ; la cinquième conserve les résultats Windows reçus et les contrôles de rapprochement. Le fichier Windows `mode3/angular_predictions.npz` possède le SHA-256 `0e6ffe660a7a8b3fbffdfe714c6d48471b0c66b8a226557fe6c97d5cb7d51d8e`.

## Compilation du PDF

Depuis `manuscript/`, avec une installation LaTeX comprenant les paquets usuels :

```sh
pdflatex -interaction=nonstopmode -halt-on-error note_pi4.tex
pdflatex -interaction=nonstopmode -halt-on-error note_pi4.tex
```

Le PDF diffusé possède le SHA-256 `f3ce473a4709f87503a7e115a26eb7d8e710d51a6423a4013296f3124ffbd9f5` et fait 537 768 octets. Ses huit pages ont été examinées individuellement à 200 dpi. Toute modification invaliderait ce contrôle et imposerait une nouvelle compilation et un nouvel examen.

## Licences et relecture

Les nouveaux textes, figures et données de Patrick Royer sont sous **CC BY 4.0** ; son code original est sous **MIT**. Voir le fichier [`../../LICENSE`](../../LICENSE) pour la portée par type de contenu et la priorité des avis propres aux tiers.

Une relecture spécialisée est sollicitée sur l’antériorité du certificat, l’interprétation des effets à taille finie et la construction d’un état réduit. Aucune expertise humaine externe n’est déclarée acquise.
