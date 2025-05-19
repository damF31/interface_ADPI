# Contexte global

Je développe une interface graphique Python avec Tkinter pour générer/éditer des fichiers d'entrée FORTRAN pour des calculs thermodynamiques d'alliages à défauts (ponctuels et complexes).  
Le projet est organisé en plusieurs fichiers (un par onglet, plus la logique principale et le multilingue).

**Le dépôt GitHub (privé/public) se trouve ici** :  
[lien_vers_ton_repo_github]  ← (remplace par l’URL de ton repo)

## Organisation des fichiers (cf. README.md pour plus de détails)

- `fortran_input_editor.py` : fichier principal, création de la fenêtre, gestionnaire des onglets, callbacks globaux, gestion de la langue, ajout/retrait dynamique des onglets.
- `tab1_general_parameters.py` : Onglet 1 – saisie des paramètres globaux (nombres de types chimiques, sous-réseaux, options défauts complexes/chargés, etc.), avec génération dynamique des champs selon les valeurs.
- `tab2_species_parameters.py` : Onglet 2 – gestion des propriétés des espèces chimiques (dynamique selon le nombre saisi en onglet 1).
- `tab3_sublattice_parameters.py` : Onglet 3 – gestion des sous-réseaux (dynamique selon le nombre saisi en onglet 1).
- `tab4_defauts_charges.py` : Onglet 4 – paramètres pour défauts ponctuels chargés (dynamique, affiché uniquement si option cochée en onglet 1).
- `tab5_complex_defects.py` : Onglet 5 – paramètres pour défauts complexes (dynamique, affiché uniquement si option cochée en onglet 1).
- `languages.py` : Dictionnaire multilingue (FR/EN) pour tous les labels.
- `aide_GIT.md` : Guide pour la gestion du dépôt avec Git/GitHub.
- (autres modules utilitaires : sauvegarde/chargement, etc.)

## Fonctionnalités clefs

- Interface multi-onglets, avec champs dynamiques selon les paramètres saisis
- Multi-langue (FR/EN) : tous les labels sont traduits dynamiquement, changement à chaud
- Onglets 4 et 5 ajoutés/retirés dynamiquement selon les cases cochées en onglet 1
- Sauvegarde/chargement des paramètres
- Vérification de cohérence minimale

---

# Ce que j’attends de toi (IA) pour cette session

- **Tu as le droit (et l’obligation) d’explorer tous les fichiers du dépôt** pour mieux comprendre le contexte.
- Prends en compte toute la structure du projet avant de proposer des modifications.
- Je veux : [décris ton besoin : debug, ajout de fonctionnalité, refactoring, documentation, autre…]
  - [exemple : “Je veux améliorer la gestion dynamique de l’onglet 5”, ou “Peux-tu trouver pourquoi le changement de langue fait planter les onglets dynamiques ?”]
- Si tu as besoin d’informations sur le fonctionnement métier (physique, chimie, modèles), demande-moi !

---

# Pour bien démarrer

Si tu as besoin d’une présentation rapide du projet, lis le fichier `README.md` à la racine du dépôt.

Merci de me guider étape par étape si tu proposes des modifications importantes, et d’indiquer précisément les fichiers/sections concernés.

