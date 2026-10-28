# Plan court de réplication indépendante — v0.2

**Base : avis DYN-VALIDEX v0.4 du 27 octobre 2026. Statut : PRÉPARÉ; PILOTE HISTORIQUE EXÉCUTÉ; COLLECTE NOUVELLE NON EXÉCUTÉE.** Cette passe ne produit ni nouvelle translation, ni graine, ni peeling scientifique nouveau, ni analyse confirmatoire. Elle a seulement rejoué les deux coordonnées historiques préspécifiées pour mesurer le coût et qualifier la seconde implémentation. Elle ne modifie pas le dossier v0.3.

## Question et décision principale

La collecte future comporte **64 lignes nouvelles complètes**. Une ligne est une translation entière avec toutes ses cellules appariées : `R∈{256,512}`, orientations `θ∈{0, atan2(1,8), π/4}` et `q_hist∈{0,95;0,90;0,85}`. Les cellules ne sont pas des réplications indépendantes.

Pour chaque cellule, avec `d=L_(1/3)-L_(0,32)`, on soustrait la sensibilité empirique

`B = 2ε_interp|a_ref,1/3-a_ref,0,32|`

`+2|a_peel-a_ref,1/3|ε_ref,1/3+ε_ref,1/3²`

`+2|a_peel-a_ref,0,32|ε_ref,0,32+ε_ref,0,32²`

`+2ε_interp(ε_ref,1/3+ε_ref,0,32)`.

Le primaire est fixé à `R=512` :

`W_j=(1/9)Σ_(3 angles × 3 q)[d_(j,θ,q)-B_(j,θ,q)]`.

La règle est `LCB = mean(W)-t_(0,975;63) sd(W)/√64 > 0`. Cette borne de Student est **nominale et de couverture approchée**; `B` est une sensibilité empirique, pas une enveloppe certifiée. Si les 64 lignes ne sont pas complètes, si les contrôles géométriques échouent ou si `LCB≤0`, aucune préférence confirmatoire n'est annoncée (`RÉPLICATION_INCOMPLÈTE` ou `NON_DISCRIMINE`). Une issue positive ne vaut que pour cette moyenne finie; elle ne dit pas que l'effet est uniforme en angle ou en contraction.

`β=0,32` reste un comparateur exploratoire fini, pas une loi, un exposant exact ou une réfutation du flot affine asymptotique.

## Secondaires préspécifiés

Ils seront tous rapportés, même défavorables, sans pouvoir sauver un primaire négatif :

1. contraste apparié `D_512-D_256`;
2. contrastes par orientation;
3. contrastes par `q`;
4. cellule dédiée `θ=π/4, q=0,95`, où les archives favorisaient `1/3`;
5. pertes complètes des sept valeurs `β∈{0,30;0,31;0,32;1/3;0,34;0,35;0,36}`, par `R`, angle et `q`.

Aucun optimum entre points, aucune phase, vitesse ou valeur de `β` ne sera ajusté après collecte. Le repère reste celui de la forme : `c_m=(2π)⁻¹∫h(φ)e^(-imφ)dφ`, `a_m=c_m/c₀`, `a_m^forme=a_m^global e^(imθ)`.

## Oracle qualifié, à geler avant collecte

Le domaine est fixé à `q_abs∈[0,75;1]`, avec

`q_abs=q_hist P_polygone_initial/(2πR)`.

Ce domaine est prudent : `h₀≥0,925` implique `B(0,925)⊂K₀`; pour toute translation, l'enveloppe de ce disque sur la grille de maille `1/R` contient `B(0,925-√2/R)`. La monotonie du périmètre donne, pour `R≥256`, `q_abs≥0,85(0,925-√2/256)>0,78` et `q_abs≤0,95`.

La qualification préparatoire a été exécutée sur 1 025 nœuds de `[0,75;1]`, sous `threadpool_limits(1)`, avec formulations, solveurs, versions, paramètres et tolérances consignés. Elle passe les contrôles de rayon de courbure, formulation distincte, raffinements séparés espace/temps, positivité de `ρ`, disques pour chaque `β` et ellipse pour `β=1/3`. Les erreurs sont des maxima empiriques sur cette grille finie, non des enveloppes de domaine ni des intervalles certifiés. Avant toute collecte, ces fichiers et leur environnement doivent encore être gelés dans le manifeste d'état.

Les futurs `q_abs` ne sont pas connus avant les enveloppes initiales. Leur évaluation mécanique par l'oracle déjà gelé est autorisée et journalisée; elle n'autorise aucun réajustement. Un échec de qualification avant collecte bloque le lancement; une sensibilité trop grande après collecte donne `NON_DISCRIMINE`.

## Aléa et exclusion exacte

Le dénominateur reste `D=2²⁰`. Les **48 couples entiers historiques** — 32 de P05 et 16 de P06–P19, intersection vide — sont repris exactement dans le JSON compagnon. P05 est un autre échantillon de la même loi dyadique que P06–P19, pas une loi différente. Les appariements utilisent `(tx,ty,D,provenance)`, jamais le numéro de ligne seul.

Une seule graine maîtresse sera fournie et consignée **au lancement futur seulement**, après gel du protocole, des programmes, de l'oracle et du manifeste, et avant toute nouvelle évaluation. Son format est une chaîne hexadécimale de 256 bits, exactement 64 caractères `[0-9a-f]`, interprétée comme un entier non signé en ordre big-endian. Aucune balise publique n'est exigée et aucune valeur n'est choisie ici.

Le générateur est `numpy.random.PCG64`, avec la version exacte de NumPy figée dans le manifeste d'environnement. À partir de l'unique graine maîtresse, le lanceur appelle séquentiellement `integers(0, 2**20, size=2, dtype=numpy.uint64)`. Chaque paire est rejetée si elle appartient à la liste noire historique ou a déjà été acceptée; les appels continuent jusqu'à 64 acceptations. L'ordre d'acceptation constitue l'ordre figé des lignes et n'est jamais trié ni permuté. Les flux auxiliaires sont initialisés uniquement par `digest=SHA-256(bytes.fromhex(master_seed_hex) || 0x00 || domain_ascii || 0x00 || bytes.fromhex(frozen_manifest_sha256))`, puis `PCG64(int.from_bytes(digest,"big"))`; les étiquettes ASCII de domaine sont distinctes et enregistrées dans le manifeste. Ces valeurs sont des dérivations déterministes de la graine maîtresse, pas de nouvelles graines indépendantes. Aucune graine ni liste n'est générée ici.

## Pilote historique exécuté et budget complet

Le pilote a rejoué seulement les coordonnées historiques `(799351,488675)/D` et `(60055,642819)/D`, sur deux rayons et trois orientations. Ces résultats sont exclus de l'analyse future des 64 lignes. Aucun tirage, aucune liste nouvelle et aucune analyse confirmatoire n'ont été exécutés.

Il a mesuré séparément, en temps mur, CPU, mémoire maximale et stockage :

- le moteur principal sur une ligne complète;
- le second moteur dans son régime réel : toutes les couches jusqu'à `q=0,85`, aux deux rayons et trois angles;
- les références continues et leurs qualifications.

Le budget livré totalise `64×coût_moteur_principal + 8×coût_second_moteur + coût_références`, et l'analogue pour le stockage. Sur les deux lignes historiques, le maximum mesuré est de 12,9446 s par ligne principale et 47,3270 s par ligne QuickHull complète. La base future séquentielle vaut 1 208,051 s; la valeur haute ×1,5, qualification incluse, vaut 1 839,868 s; le stockage haut ×1,25 vaut 259 226 302,5 octets. Ces chiffres n'emportent aucun GO. Les huit lignes du second moteur restent fixées par leur position dans la future liste : `0,8,16,24,32,40,48,56`.

La marche de Jarvis exhaustive a été projetée à au plus 101 475,60 s par ligne complète et n'a pas été lancée. Le contrôle retenu est un QuickHull entier complet avec boucle de couches et retrait distincts; sa limite déclarée est le partage des intervalles initiaux qualifiés. Le plan courant reste `n=64`; si ce coût est inacceptable, une nouvelle version prospective chiffrée doit être décidée avant tout tirage. L'effectif ne sera jamais choisi après lecture de nouvelles pertes.

## Intégrité et arrêt

Les deux moteurs retirent uniquement les sommets extrêmes et conservent les points alignés non sommets. Les contrôles portent sur l'appartenance initiale, chaque couche utilisée, les enveloppes adjacentes, l'imbrication, les poids de Minkowski et les coefficients issus des mêmes primitives. Une divergence de sommets bloque l'analyse et ne justifie jamais la suppression d'une ligne.

Tout correctif géométrique est aveugle aux contrastes : mêmes coordonnées, aucune exclusion, nouvelle empreinte du code et recalcul de tous les résultats dépendants avant analyse. Pas d'arrêt anticipé pour efficacité ou futilité; reprise seulement de la même coordonnée en cas de panne technique. Une campagne incomplète n'est pas analysée comme si `n` avait été réduit.

## Deux corrections documentaires à reporter

P05 est un autre échantillon de la même loi dyadique que P06–P19, et non une loi différente. Seules trois trajectoires ont neuf cellules éligibles après exclusion; les treize autres du groupe historique concerné en ont dix.
