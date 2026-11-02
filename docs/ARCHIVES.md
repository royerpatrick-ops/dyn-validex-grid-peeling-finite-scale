# Archives de production et de provenance

| Fichier | Taille (octets) | SHA-256 | Rôle |
|---|---:|---|---|
| `_prepare_replication_v0_4_20261027.zip` | 210150972 | `f225c317be686798471e52a2f34e2bbc2fdf22406654702aeebf197797764dc3` | Production après collecte, 64 sauvegardes géométriques |
| `_execute_replication_v0_5_20261027.zip` | 31186 | `a883fcd2776a4f3b70e0992a5312c4a2a49cf4860dfa483f35539c50e0c32378` | Traces et compte rendu d'exécution |
| `DYN-Preparation_Replication_Independante_v0_4.zip` | 7290136 | `22c2aa3efe02113feb452c48d9a5684961bfdf36a49367d640d0be2df8a28ad9` | État antérieur à la collecte |

Les trois ZIP ont passé un contrôle CRC intégral. La grande archive contient exactement les 64 fichiers principaux `row_00_main.npz` à `row_63_main.npz`, sans manque ni doublon.

La première archive n'est pas stockée dans l'historique Git. Les trois archives sont destinées aux actifs de release et au supplément Zenodo. Ne pas confondre la préparation antérieure à la collecte avec l'archive de production après campagne.

## Statut historique et limites de portabilité

Ces ZIP sont des artefacts scientifiques gelés, conservés octet pour octet ; ils ne sont pas des paquets publics nettoyés. Ils contiennent des chemins absolus historiques sous `C:\Users\dynagenese\...`. L'archive de production contient en outre six fichiers `.pyc` dont le nom de source compilé est absolu, un ancien prompt de travail, un journal de commandes et un PDF interne antérieur. Les mentions « aucune publication » ou équivalentes décrivent uniquement l'ancienne passe v0.4 et ne remplacent pas l'autorisation de publication v0.7.

L'ancien PDF relatif à la π-morphodynamique n'est pas le manuscrit de la présente note. Les fichiers `authorization_draw_v05.json` et `authorization_run_v05.json` sont des garde-fous locaux de campagne, non des justificatifs d'autorisation GitHub ou Zenodo. Aucune licence n'est incorporée dans les trois ZIP ; les licences applicables sont fournies à l'extérieur, dans le présent dépôt. Les scripts et manifestes historiques à chemins absolus ne doivent pas être exécutés aveuglément : utiliser les données compactes et les adaptateurs documentés dans le dépôt.

Ces traces ont été signalées plutôt que réécrites afin de préserver leurs empreintes d'origine. Aucun secret d'authentification, mot de passe, jeton, clé privée, adresse électronique réelle ou article tiers n'a été détecté lors de l'audit avant diffusion.

## Complément π/4 — v0.2

La release v0.2 et la version Zenodo correspondante distribuent cinq ZIP supplémentaires, sans les ajouter à l'historique Git :

| Fichier | Taille (octets) | SHA-256 | Rôle |
|---|---:|---|---|
| `DYN-Controle_Symetries_Pi4_v0_1.zip` | 18945079 | `028b63a97a8b01187a04be840a4ca778aa33b497114ed7fd6b230140da520242` | 18 trajectoires de contrôle supplémentaires et contrôles de symétrie |
| `DYN-Localisation_Angulaire_Mode3_Pi4_v0_1.zip` | 19507630 | `ddb3429d4fea4624b14d1a9cfacebcfcb1bab64b2c918a48dc78146ebca085f1` | localisation signée du mode 3 et compensations |
| `DYN-Deplacements_Repere_Fixe_Pi4_v0_1.zip` | 66232366 | `0898cb5a9cc6ceb52929bbaeec7db24ec72320acc9be0ffdfd2869635c0f566f` | déplacements à repère fixe et référence affine sauvegardée |
| `DYN-Sauts_Mode3_Pi4_Prediction_Conditionnelle_v0_1.zip` | 25709536 | `d40985bfb59fce13441a049eb2a882cf8ce82deb9be5be9be70238a8b347099e` | certificats de sauts et prévision conditionnelle à un pas |
| `DYN-Pi4_Resultats_Windows_Verifies_20261102.zip` | 24893952 | `ea1364e07c00110ad898e924bc7a82eb5a801c3983c996a5a2bfb3deedbc6956` | résultats Windows reçus et rapprochements |

Les cinq CRC et les 181 charges utiles décrites par les manifestes internes ont été vérifiés. Les archives restent inchangées. La sélection versionnée sous `supplements/pi4/reproducibility/` ne remplace pas ces originaux : elle regroupe seulement les scripts portables et les entrées compactes utiles aux replays documentés.
