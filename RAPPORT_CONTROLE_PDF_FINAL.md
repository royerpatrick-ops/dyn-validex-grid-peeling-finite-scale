# Rapport de contrôle PDF final

Date du contrôle : 28 octobre 2026. Le protocole `PROTOCOLE_CONTROLE_PDF.md` de la transmission v0.7 a été appliqué aux deux PDF destinés au public. Le verdict porte sur les octets identifiés ci-dessous ; une modification ultérieure l'invalide.

## Manuscrit autonome v0.1

- Fichier : `output/pdf/DYN-Note_Experimentale_Grid_Peeling_v0_1.pdf`
- SHA-256 : `6a66e6a69e4fcf9c8ff21ea437c7d1ad0d2cf3c2cf2ca716e1574032ba197fbd`
- Taille : 594 872 octets ; 10 pages A4 ; PDF 1.5.
- Moteur : MiKTeX-pdfTeX 4.23 (MiKTeX 25.12 ; producteur pdfTeX-1.40.28) et MiKTeX-BibTeX 4.2.
- Polices : familles Latin Modern texte, mono et mathématiques, toutes incorporées, sous-ensemblées et munies d'un mappage Unicode.
- Journal final : aucun débordement, sous-remplissage, avertissement LaTeX/de paquet, citation ou référence non définie.
- Rendu 200 dpi : pages 1 à 10 rendues et examinées individuellement.
- Rendu 300 dpi : pages 1 à 10 rendues ; pages 4, 5, 7, 9 et 10 examinées en détail pour les formules, tableaux, légendes, références et empreintes.
- Extraction : DOI Zenodo et DOI des deux références retrouvés exactement ; aucun caractère de remplacement.

Les contrôles intermédiaires avaient révélé des débordements d'empreintes SHA-256, une hauteur d'en-tête insuffisante, des définitions et une liste de dépendances trop larges, une collision entre le nom de la macro DOI et la bibliographie, ainsi qu'un placement de figure inadéquat. Ces cinq catégories ont été corrigées avant la compilation finale. La version portant l'empreinte ci-dessus ne présente aucune collision, superposition, découpe, formule débordante, glyphe manquant, caractère doublé, tableau illisible ou référence orpheline.

Verdict de mise en page : **PASS**.

## Fiche de synthèse historique v0.6

- Fichier : `evidence/v0.6/DYN-Fiche_Replication_Prospective_64_v0_6.pdf`
- SHA-256 : `6249de781857c9a97bb181faee9f21d52d66191d47025a4a775c0748de871084`
- Taille : 142 621 octets ; 3 pages A4 ; PDF 1.4.
- Producteur : ReportLab PDF Library, sans version inscrite dans le fichier.
- Polices : DejaVu Sans et DejaVu Sans Bold incorporées/sous-ensemblées avec mappage Unicode ; Helvetica Type 1 non incorporée et sans table ToUnicode.
- Rendu 200 et 300 dpi : pages 1, 2 et 3 rendues puis examinées individuellement aux deux résolutions.
- Extraction : les glyphes π, ρ, φ, β, θ, ×, ₀ et ₃ sont retrouvés ; aucun caractère de remplacement.

`pdffonts` a signalé l'absence de polices d'affichage locales pour plusieurs alias. Cette pièce historique a néanmoins été rendue correctement par Poppler : aucune collision, superposition, découpe, duplication, perte de glyphe ni illisibilité n'a été observée. Elle est conservée octet pour octet ; aucun correctif n'y a été appliqué. Helvetica restant substituable, la réserve de portabilité est explicite.

Verdict de mise en page : **PASS avec réserve technique sur l'incorporation de Helvetica**.

## Distinction scientifique

Ces verdicts évaluent la lisibilité et l'intégrité de mise en page. Ils ne lèvent pas les limites scientifiques : borne Student nominale et approchée, budgets numériques empiriques, contrôle de 30 863 couches limité à huit lignes et portée expérimentale restreinte. La fiche v0.6 reste une synthèse contrôlée, non un substitut au manuscrit autonome.

Verdict global : **PASS, avec la réserve de police documentée pour la fiche historique v0.6**.
