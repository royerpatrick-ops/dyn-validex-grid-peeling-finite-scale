# Provenance des dépendances gelées retrouvées

| Origine | Chemin public relatif | SHA-256 |
|---|---|---|
| P06 | `src/recovered_dependencies/P06/support_shape_engine.py` | `f1769af06d914f5ff078a3c22c5c259b8777656d5c3596321362766a6e51a191` |
| P07 | `src/recovered_dependencies/P07/baseline/polygon_fourier.py` | `ec5bcec9b6e5c52bc136f39ee9c27a82a2fbc7c03b2493b24de321c3b4446245` |
| P07 | `src/recovered_dependencies/P07/targeted_supports.npz` | `5efc10d213e334b0a09d8686b4f7aa3f4f6950d6735a96cbba6999165517530e` |
| P08 | `src/recovered_dependencies/P08/extension_supports.npz` | `36de3f31f8fdd3c2f02e680a18433f272623c9abfbbf6aeac8452307d69e3ea4` |

Les empreintes correspondent exactement à `SOURCES_ET_DEPENDANCES.json`. Les deux NPZ passent leur CRC interne. Les protocoles et manifestes P06/P07/P08 chaînent ces fichiers comme entrées gelées.

Aucun fichier LICENSE, avis SPDX ou en-tête de licence antérieur n'a été retrouvé dans ces trois paquets. La publication ne prétend donc pas avoir découvert une licence préexistante ; la portée de diffusion décidée par l'auteur est documentée dans `LICENSE`.
