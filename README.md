# ==============================================================================
# interface_ADPI
# ==============================================================================

# Interface graphique ADPI – Générateur de fichiers d'entrée FORTRAN pour la thermo

Ce projet est une interface graphique en Python (Tkinter) pour créer, éditer et sauvegarder des fichiers d'entrée à destination de codes FORTRAN de calculs thermodynamiques sur les alliages à défauts.  
L’objectif est de faciliter la saisie, la modification et la validation de ces fichiers complexes via une interface claire, structurée en onglets thématiques.

---

## Fonctionnalités principales

- **Interface multi-onglets** pour organiser la saisie par grands thèmes
- **Champs dynamiques** : le nombre de champs dépend des paramètres saisis (ex : nombre de sous-réseaux, espèces…)
- **Multi-langue (français/anglais)** : tous les libellés sont traduits dynamiquement
- **Onglets ajoutés/retirés dynamiquement** selon les options (défauts chargés, complexes…)
- **Sauvegarde et chargement** des paramètres (format texte ou JSON)
- **Vérification de cohérence** basique sur la saisie

---

## Structure des onglets

### **Onglet 1 – Paramètres généraux (‘Tab1GeneralParameters’)**
- Commentaire
- Nom de base des fichiers de sortie
- Nombre de types chimiques (total & intrinsèque)
- Nombre de sous-réseaux
- Présence d’interstitiels et d’interstitiels intrinsèques
- Nombre de sous-réseaux occupés (par espèce intrinsèque) – dynamique
- Nombre de sites par maille (par sous-réseau) – dynamique
- Activation défauts complexes (active l’onglet 5)
- Activation défauts ponctuels chargés (active l’onglet 4)
- Température, pression
- Énergies de référence (par type chimique) – dynamique
- Mode thermodynamique (par atome / par maille)
- Type d’énergie libre (totale / formation)
- Type de calcul (muVT, NPT, NPT=0)
- etc.

### **Onglet 2 – Paramètres espèces chimiques (‘Tab2SpeciesParameters’)**
- Saisie détaillée des propriétés de chaque espèce chimique intrinsèque ou interstitielle
- Nom ou symbole chimique
- Masse molaire
- Numéro d’ordre
- Statut (ex : “actif”, “interstitiel”, “solvant”, etc.)
- (Éventuellement) Couleur ou code de représentation graphique
- **Champs dynamiques :** nombre d’espèces = paramètre défini dans l’onglet 1

### **Onglet 3 – Paramètres sous-réseaux (‘Tab3SublatticeParameters’)**
- Saisie des propriétés pour chaque sous-réseau
- Nom ou numéro du sous-réseau
- Type (ex : métallique, interstitiel…)
- Espèces autorisées sur chaque sous-réseau (liste déroulante ou cases à cocher)
- Nombre de sites par maille (recopie/édition ?)
- (Optionnel) Paramètres structuraux spécifiques
- **Champs dynamiques :** nombre de sous-réseaux = paramètre défini dans l’onglet 1

### **Onglet 4 – Défauts ponctuels chargés (‘Tab4DefautsCharges’)**
- Apparaît seulement si l’option “défauts ponctuels chargés” est cochée dans l’onglet 1
- Paramètres relatifs aux défauts chargés
- Exemple : niveau de Fermi, corrections, état de charge possible, etc.
- Adaptable selon les besoins du calcul

### **Onglet 5 – Défauts complexes (‘Tab5ComplexDefects’)**
- Apparaît seulement si l’option “défauts complexes” est cochée dans l’onglet 1
- Nombre maximum de sites des complexes
- Nombre de types de défauts complexes
- Pour chaque complexe :
    - Numéro, multiplicité d’orientation, sous-réseau associé, nombre de sites
    - Numéros de sous-réseaux occupés (liste dynamique)
    - Types chimiques sur chaque site (liste dynamique)
    - Énergie, volume de supercellule
    - (Commentaire ou label facultatif)

---

## Structure des fichiers du projet

- **`fortran_input_editor.py`** : Fichier principal, lance la fenêtre, instancie et gère tous les onglets, gestion des changements de langue, callbacks globaux.
- **`tab1_general_parameters.py`** : Onglet 1, logique des paramètres généraux, callbacks pour l’ajout/retrait dynamique d’onglets.
- **`tab2_species_parameters.py`** : Onglet 2, gestion des champs dynamiques pour chaque espèce chimique.
- **`tab3_sublattice_parameters.py`** : Onglet 3, gestion des sous-réseaux, champs dynamiques pour chaque sous-réseau.
- **`tab4_defauts_charges.py`** : Onglet 4, spécifique aux défauts ponctuels chargés (optionnel, dynamique).
- **`tab5_complex_defects.py`** : Onglet 5, spécifique aux défauts complexes (optionnel, dynamique).
- **`languages.py`** : Dictionnaire de traduction FR/EN pour tous les labels de l’interface.
- **`aide_GIT.md`** : Guide pour la gestion du projet avec Git/GitHub.
- (Autres : utilitaires, modules de sauvegarde/chargement…)

---

## Avancement / TODO

- [x] Onglets 1 à 3 : champs dynamiques et navigation fonctionnelle
- [x] Dynamique d’ajout/retrait des onglets 4 et 5 selon les options
- [x] Gestion multi-langue
- [x] Sauvegarde/chargement basique OK
- [x] Gestion dynamique des défauts complexes (onglet 5)
- [ ] Compléter et valider les comportements avancés de l’onglet 4 (défauts chargés)
- [ ] Améliorer la robustesse à la reconstruction de l’interface (changement de langue, etc.)
- [ ] Ajout de tests, doc utilisateur avancée, gestion d’aide contextuelle

---

## Lancement

```bash
python main_v2.py
```

---

## Dépendances

- Python ≥ 3.8
- Tkinter (livré de base)
- (optionnel : numpy, pandas…)

---

## Pour contribuer / demander de l’aide

Voir le fichier `aide_GIT.md` pour la gestion du dépôt.  
Si tu demandes de l’aide à une IA, donne le lien du dépôt et précise le ou les fichiers concernés, ainsi que la partie de l’interface ou la logique à modifier.

---

