# π-morphodynamique : réplication prospective

64 nouvelles translations - Fiche DYN-VALIDEX v0.6 - 27 octobre 2026
Patrick Royer, chercheur indépendant non affilié - assistance IA

## Résultat principal : le critère fixé avant collecte est satisfait

La campagne unique est complète et son résultat primaire est reproduit indépendamment : à R=512, le modèle continu de vitesse normale κ^0,32 a une perte moyenne sur le coefficient de forme a₃ inférieure à celle du flot affine κ^(1/3), après soustraction de la sensibilité numérique empirique prescrite. La décision du programme gelé est FINITE_R512_PREFERENCE_0_32.

Ce résultat concerne une comparaison moyenne à taille finie. Il ne détermine pas un exposant exact et ne démontre ni uniformité angulaire, ni convergence, ni réfutation d'une limite affine. Les résultats secondaires retrouvent précisément l'hétérogénéité observée pendant l'exploration.

## Cadre et plan conservés

La forme initiale a pour fonction d'appui h₀(φ)=1-(3/40)cos(3φ), avec rayon de courbure ρ₀=1+(3/5)cos(3φ)>0. Les deux rayons sont R=256 et 512 ; les orientations sont 0, atan(1/8) et π/4 ; les contractions historiques sont q=0,95, 0,90 et 0,85. Les références sont comparées au périmètre absolu propre à chaque cible, sans ajustement libre d'horloge ou de phase.

Le plan et le code v0.4 ont été gelés avant la graine et la collecte. Les 64 nouvelles translations sont distinctes, tirées sans remise sur la grille dyadique de dénominateur 2^20, hors des 48 paires historiques. Les mêmes translations sont appariées entre rayons, angles, contractions et références. Les anciennes 16 translations exploratoires ne sont pas ajoutées au test confirmatoire.

On note D=L_(1/3)-L_(0,32), où L est la perte |a₃_peel-a₃_ref|² et a₃=c₃/c₀. La sensibilité B comprend l'interpolation, les erreurs empiriques des références et leur terme mixte. Le score principal par translation est W_j=moyenne_(3 angles × 3 q)(D-B), à R=512. Le test fixé est une borne Student nominale unilatérale à 97,5 %, sur les 64 scores W.

| Quantité primaire recomposée | Valeur |
|---|---|
| Nombre de translations nouvelles | 64 |
| Moyenne W, après soustraction de B | 1,1220297693×10^-7 |
| Écart-type entre translations | 8,1552119603×10^-9 |
| Quantile t_(0,975;63) | 1,9983405425 |
| Borne inférieure nominale | 1,1016586560×10^-7 > 0 |
| Scores W positifs | 64 sur 64 |
| Bornes positives après retrait d'une translation | 64 sur 64 |

La couverture Student est nominale et approchée. Les sensibilités B sont empiriques, non certifiées par intervalles. Le succès est celui du critère prévu, dans ce cadre numérique explicite.

<!-- PAGEBREAK -->

## Résultats secondaires : déplacement avec le rayon et hétérogénéité

Les pertes moyennes suivantes sont calculées sur les neuf cellules de chaque translation puis sur les 64 translations. Les 576 valeurs par rayon ne constituent pas 576 réplications indépendantes.

| Rayon | Perte β=1/3 | Perte β=0,32 |
|---|---|---|
| 256 | 6,53293×10^-8 | 2,06660×10^-7 |
| 512 | 2,37774×10^-7 | 7,94367×10^-8 |

Parmi les sept β prévus, 1/3 a la plus petite perte agrégée à R=256 et 0,32 à R=512. Le contraste brut D512-D256 vaut en moyenne 2,99667823×10^-7 et est positif pour les 64 translations. Ce sont des secondaires préspécifiés ; ils ne remplacent pas le primaire et n'identifient aucun optimum entre les points de grille.

À R=512, huit des neuf moyennes brutes D sont positives, mais six seulement gardent une marge D-B positive. La cellule θ=π/4, q=0,95 favorise 1/3 pour les 64 nouvelles translations, avec D moyen=-4,34163977×10^-8. La préférence moyenne pour 0,32 ne vaut donc pas uniformément sur la grille observée.

![Contrastes par orientation et contraction à R=512](contrastes_R512.png)

Figure : moyennes sur les 64 translations, en unités de 10^-7. À gauche, contraste brut D ; à droite, marge conservatrice D-B. Une marge négative à droite ne démontre pas, à elle seule, une préférence pour 1/3 ; elle indique l'absence de marge favorable à 0,32 sous la sensibilité soustraite.

La perte affine agrégée augmente d'un facteur 3,63963 de R=256 à 512. L'expérience ne montre donc pas une amélioration monotone avec le raffinement sur ce couple de rayons. Cette observation finie ne décide pas une limite asymptotique.

<!-- PAGEBREAK -->

## Vérifications indépendantes réalisées sur les pièces reçues

Les deux archives passent le CRC ; les 24 fichiers gelés livrés sont inchangés. Les traces relient le même manifeste au tirage, à la campagne et à l'analyse. Le rejeu de la graine sauvegardée reproduit les 64 paires et leur ordre, sans doublon ni intersection avec les 48 paires historiques.

Le calcul autonome, sans importer ni réexécuter l'analyseur gelé, retrouve les 8 064 cellules, 1 152 groupes et 64 scores. Écart maximal des statistiques primaires : 3,31×10^-24. Les références et sensibilités concordent avec les NPZ. Les 1 152 coefficients interpolés a₃ ont été reconstruits depuis les sommets des 64 NPZ ; écart maximal au CSV : 7,62×10^-17.

La géométrie sauvegardée des huit lignes prescrites a été contrôlée sur 30 863 couches : intervalles entiers, retraits des seuls sommets, comptes, convexité stricte, imbrication, encadrements et digests concordent. L'intégration de Fourier indépendante donne des écarts maximaux de 3,26×10^-16 aux extrémités et 2,98×10^-16 après interpolation. Aucun nouveau hull ni épluchage n'a été calculé lors de ces contrôles.

Le contrôle géométrique complet porte sur huit lignes, pas 64. L'appartenance initiale au corps lisse n'a pas été requalifiée ici ; QuickHull partage certaines primitives avec le moteur principal et ses polygones ne sont pas sauvegardés séparément. Quatre dépendances originales gelées restent hors des copies livrées.

## Coût réalisé et portée

RUN a pris 1 216,40 secondes, soit 20 minutes 16,4 secondes. La campagne et l'analyse occupent ensemble environ 208,4 Mo, sous les prévisions. Le pic mémoire de production n'a pas été enregistré. État final : COMPLETE ; aucune seconde campagne ni publication automatique n'est rapportée.

Le contraste fini préspecifié est répliqué sous les conventions annoncées. Les nouvelles translations séparent exploration et test prospectif ; elles ne fournissent ni expertise humaine indépendante, ni certification par intervalles.

Le lien épluchage-flot affine est connu. L'objet documenté ici est la comparaison finie, appariée, à périmètre égal ; la nouveauté bibliographique n'est pas évaluée.

## Décision de suite

Le jalon prévu est terminé. La prochaine étape utile est une note expérimentale pour relecture spécialisée, conservant le primaire positif et les résultats angulaires défavorables. Exploration et 64 nouvelles translations doivent rester séparées. Une autre étude de rayons ou de formes nécessiterait une nouvelle question précise.

La preuve statique H_rot, l'hypothèse à orientation fixe et toute relation générale avec « π efficace » demeurent des objets distincts. Cette réplication n'en modifie pas le statut et n'introduit aucune nouvelle valeur de π.

## Sources et traces

Sources reçues : _execute_replication_v0_5_20261027.zip et _prepare_replication_v0_4_20261027.zip. Empreintes, scripts et bilans autonomes sont joints séparément.

Contexte primaire déjà consulté : Eppstein-Har-Peled-Nivasch, https://arxiv.org/abs/1710.03960 ; Rote-Rüber-Saghafian, https://doi.org/10.4230/LIPIcs.SoCG.2024.76.
