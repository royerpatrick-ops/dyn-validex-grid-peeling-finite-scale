# Replays portables du complément π/4

Ce dossier est une sélection byte-identique de scripts et d’entrées compactes provenant de `DYN-Sauts_Mode3_Pi4_Prediction_Conditionnelle_v0_1.zip` (SHA-256 `d40985bfb59fce13441a049eb2a882cf8ce82deb9be5be9be70238a8b347099e`). L’archive originale complète reste l’autorité de provenance et est distribuée avec la release v0.2 et Zenodo.

Python, NumPy et SciPy sont requis ; les versions qualifiées sont consignées dans `requirements.txt`.

## Certificats d’amplitude

```powershell
python run_amplitude_checks.py --output amplitude_replay
```

Résultat attendu : `PASS_AMPLITUDE_CERTIFICATES`, avec 258 077 cas certifiés complets, 4 584 censurés, 703 abstentions informatives et aucun échec certifié complet. Le replay lit les sections courantes extraites et les trajectoires de support sauvegardées ; il ne relance pas l’extraction géométrique.

Les résultats du replay Windows pour ce certificat d’amplitude n’ont pas été reçus. La preuve et les calculs de référence fournis restent disponibles ; cette limite est documentaire et non un motif de nouvelle campagne.

## Prévision conditionnelle globale du mode 3

```powershell
python run_mode3_prediction.py --output results_replay
```

Résultat attendu : `PASS_CONDITIONAL_ENDPOINT_MODE3_PREDICTION`, 1 152 prédictions, aucune divergence géométrique ni de population et une erreur maximale `a3` de l’ordre de l’arrondi machine.

Le script écrit et scelle les prévisions avant de lire le benchmark. Il reçoit le polygone complet courant et applique un pas du peeling microscopique ; il ne prédit ni l’évolution depuis la forme initiale ni le temps d’atteinte d’une cible.

## Diagnostic à seize directions

```powershell
python run_sparse_angular_checks.py --output sparse_replay
```

Résultat attendu : `DIAGNOSTIC_NOT_FORECAST`. Les seize supports sont observés, non prédits. Le script construit un autre polygone convexe partageant ces échantillons et mesure la perte d’information spatiale ; l’échec des deux reconstructions spécifiées ne réfute pas toute fermeture possible.

## Limites d’exécution

- Utiliser un nom de sortie neuf ; les scripts refusent d’écraser un dossier existant.
- Aucun script de ce dossier ne lance de tirage, de campagne libre, de trajectoire de peeling, de solveur EDP ou d’ajustement.
- Les scripts historiques à chemins absolus ne sont pas copiés ici ; ils restent inchangés dans les archives de provenance.
- Les sorties de replay sont des reçus locaux et ne doivent pas être ajoutées au dépôt sans nouveau protocole de publication.
