# Guide rapide : Gérer ce projet avec GIT et GitHub

## 1. Créer un dépôt GitHub et initialiser le projet

### a) Sur GitHub.com

1. Clique sur “New repository” ou “Nouveau dépôt”
2. Donne un nom (ex : `interface_ADPI`)
3. Choisis “Public” ou “Private”
4. (Optionnel) Ajoute un README.md initial
5. Clique sur “Create repository”

### b) Sur ta machine

1. Place-toi dans le dossier du projet (ou clone le dépôt si tu pars de zéro) :
    ```bash
    cd /chemin/vers/ton/projet
    ```
2. Si ce n’est pas déjà un dépôt git, initialise-le :
    ```bash
    git init
    ```
3. Lien avec GitHub :
    ```bash
    git remote add origin https://github.com/ton_login/nom_du_repo.git
OU
    git remote set-url origin git@github.com:damF31/interface_ADPI.git
    ```

## 2. Ajouter et pousser tes fichiers

1. Ajoute tous tes fichiers au suivi git :
    ```bash
    git add .
    ```
2. Fais un commit avec un message :
    ```bash
    git commit -m "Premier commit : ajout de la base du projet"
    ```
3. Pousse sur GitHub :
    ```bash
    git push -u origin master
    ```
   (ou `main` selon le nom de ta branche par défaut)

## 3. Mettre à jour le dépôt

1. Vérifie les fichiers modifiés :
    ```bash
    git status
    ```
2. Ajoute ceux à pousser :
    ```bash
    git add fichier_modifié.py
    ```
3. Committe et pousse :
    ```bash
    git commit -m "Description des changements"
    git push
    ```

## 4. Récupérer les nouveautés du dépôt distant (GitHub)

```bash
git pull
```

## 5. Conseils pratiques

- **Commite souvent**, avec des messages clairs.
- Utilise `.gitignore` pour ne pas versionner les fichiers temporaires ou spécifiques à ta machine.
- En cas de conflit lors d’un `git pull`, lis bien les messages : git te guide pour résoudre les conflits.
- Pour demander de l’aide à l’IA, donne le lien du dépôt et explique ce que tu veux faire ou corriger.

## 6. Ressources utiles

- [Documentation Git officielle](https://git-scm.com/doc)
- [Guide GitHub débutant](https://docs.github.com/fr/get-started/quickstart)
- [Aide mémoire Git (fr)](https://rogerdudler.github.io/git-guide/index.fr.html)

