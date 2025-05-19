import json
from tkinter import filedialog, messagebox

# Pour l'export FORTRAN : à adapter selon ton format de fichier cible
def export_fortran(data_dict, parent_window=None):
    """Exporte les paramètres dans un fichier au format FORTRAN (texte)."""
    try:
        file_path = filedialog.asksaveasfilename(
            defaultextension=".dat",
            filetypes=[("Fichiers d'entrée FORTRAN", "*.dat"), ("Tous les fichiers", "*.*")],
            parent=parent_window
        )
        if not file_path:
            return  # Annulé

        with open(file_path, "w") as f:
            # TODO : Adapter ce bloc pour exporter chaque onglet/paramètre au bon format
            f.write("# Fichier généré par l'interface ADPI\n")
            for key, value in data_dict.items():
                f.write(f"{key} = {value}\n")

        if parent_window:
            messagebox.showinfo("Export réussi", f"Fichier exporté : {file_path}", parent=parent_window)
    except Exception as e:
        if parent_window:
            messagebox.showerror("Erreur d'export", str(e), parent=parent_window)
        else:
            print("Erreur d'export :", e)

def import_fortran(parent_window=None):
    """ImporTe un fichier FORTRAN et retourne un dictionnaire des paramètres."""
    try:
        file_path = filedialog.askopenfilename(
            filetypes=[("Fichiers d'entrée FORTRAN", "*.dat"), ("Tous les fichiers", "*.*")],
            parent=parent_window
        )
        if not file_path:
            return None  # Annulé

        data_dict = {}
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" in line:
                    key, value = line.split("=", 1)
                    data_dict[key.strip()] = value.strip()
        return data_dict
    except Exception as e:
        if parent_window:
            messagebox.showerror("Erreur d'import", str(e), parent=parent_window)
        else:
            print("Erreur d'import :", e)
        return None

def save_session(data_dict, parent_window=None):
    """Sauvegarde la session utilisateur au format JSON (état complet de l'UI)."""
    try:
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("Fichiers de session JSON", "*.json"), ("Tous les fichiers", "*.*")],
            parent=parent_window
        )
        if not file_path:
            return
        with open(file_path, "w") as f:
            json.dump(data_dict, f, indent=2)
        if parent_window:
            messagebox.showinfo("Sauvegarde réussie", f"Session enregistrée : {file_path}", parent=parent_window)
    except Exception as e:
        if parent_window:
            messagebox.showerror("Erreur de sauvegarde", str(e), parent=parent_window)
        else:
            print("Erreur de sauvegarde :", e)

def load_session(parent_window=None):
    """Charge une session utilisateur au format JSON (état complet de l'UI)."""
    try:
        file_path = filedialog.askopenfilename(
            filetypes=[("Fichiers de session JSON", "*.json"), ("Tous les fichiers", "*.*")],
            parent=parent_window
        )
        if not file_path:
            return None
        with open(file_path, "r") as f:
            data_dict = json.load(f)
        return data_dict
    except Exception as e:
        if parent_window:
            messagebox.showerror("Erreur de chargement", str(e), parent=parent_window)
        else:
            print("Erreur de chargement :", e)
        return None

# Astuce : tu peux ajouter ici d'autres fonctions utilitaires si besoin
