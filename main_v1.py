import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import json
import os

from languages import LANGUAGES
from tab1_general_parameters import Tab1GeneralParameters
from tab2_muVT_NPT import Tab2MuVTNPT
from tab3_supercell import Tab3Supercell
from file_io import export_fortran, import_fortran, save_session, load_session

class FortranInputEditor(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Fortran Input File Editor")
        self.geometry("900x700")
        self.current_language = 'fr'
        self.lang_dict = LANGUAGES[self.current_language]
        self.tabs = None
        self.create_widgets()

    def get_global_state(self):
        return {
            "chem_types_total": self.tab1.chem_types_total_var.get()
        }

    def create_widgets(self):
        if hasattr(self, "tabs") and self.tabs is not None:
            self.tabs.destroy()
        self.menubar = tk.Menu(self)
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label=self.lang_dict['open'], command=self.open_file)
        filemenu.add_command(label="Importer .dat/.inp...", command=self.open_fortran_input)
        filemenu.add_command(label=self.lang_dict['save'], command=self.save_file)
        filemenu.add_separator()
        filemenu.add_command(label=self.lang_dict['quit'], command=self.quit)
        self.menubar.add_cascade(label=self.lang_dict['file'], menu=filemenu)
        langmenu = tk.Menu(self.menubar, tearoff=0)
        langmenu.add_command(label=self.lang_dict['french'], command=lambda: self.set_language('fr'))
        langmenu.add_command(label=self.lang_dict['english'], command=lambda: self.set_language('en'))
        self.menubar.add_cascade(label=self.lang_dict['language'], menu=langmenu)
        self.config(menu=self.menubar)

        self.tabs = ttk.Notebook(self)
        self.tab1 = Tab1GeneralParameters(self.tabs, self.lang_dict, tab_manager=self)
        self.tab2 = Tab2MuVTNPT(self.tabs, self.lang_dict, self.get_global_state)
        self.tab3 = Tab3Supercell(self.tabs, self.lang_dict)
        self.tabs.add(self.tab1, text=self.lang_dict['tab1'])
        self.tabs.add(self.tab2, text=self.lang_dict['tab2'])
        self.tabs.add(self.tab3, text=self.lang_dict['tab3'])
        self.tabs.pack(expand=1, fill="both")
        self.tab4 = None  # Pour l’onglet défauts chargés
        self.tab4_exists = False

    def add_tab4(self):
       if not hasattr(self, 'tab4_exists') or not self.tab4_exists:
           from tab4_defauts_charges import Tab4DefautsCharges
           self.tab4 = Tab4DefautsCharges(self.tabs, self.lang_dict)
           self.tabs.add(self.tab4, text=self.lang_dict.get('tab4', "Options défauts chargés"))
           self.tab4_exists = True
    
    def remove_tab4(self):
       if hasattr(self, 'tab4_exists') and self.tab4_exists and hasattr(self, 'tab4'):
           idx = self.tabs.index(self.tab4)
           self.tabs.forget(idx)
           self.tab4 = None
           self.tab4_exists = False

    def add_tab5(self):
        if not hasattr(self, 'tab5_exists') or not self.tab5_exists:
            from tab5_complex_defects import Tab5ComplexDefects
            self.tab5 = Tab5ComplexDefects(self.tabs, self.lang_dict)
            self.tabs.add(self.tab5, text=self.lang_dict.get('tab5', "Défauts complexes"))
            self.tab5_exists = True
    
    def remove_tab5(self):
         # ### Correction: vérifier si self.tab5 existe et est bien un onglet du notebook
         if hasattr(self, 'tab5_exists') and self.tab5_exists and hasattr(self, 'tab5'):
            try:
                idx = self.tabs.index(self.tab5)
                self.tabs.forget(idx)
            except Exception as e:  ### <-- Ajout : gestion de l'exception si l'onglet n'existe plus
                pass
            self.tab5 = None
            self.tab5_exists = False

    def set_language(self, lang_code):
        # On sauvegarde les données pour ne pas les perdre pendant le switch
        data = self.get_all_data()
        self.current_language = lang_code
        self.lang_dict = LANGUAGES[lang_code]
        self.create_widgets()
        self.set_all_data(data)
        # ===> AJOUTE ICI le bloc de rétablissement des onglets dynamiques
        # --- Pour l'onglet 4 (défauts ponctuels chargés) ---
        if hasattr(self.tab1, 'defaut_charge_var') and self.tab1.defaut_charge_var.get():
            self.add_tab4()     ### Ajout pour restaurer l'onglet 4 si coché
        else:
            self.remove_tab4()  ### Ajout pour retirer l'onglet 4 si décoché
        # --- Pour l'onglet 5 (défauts complexes) ---
        if hasattr(self.tab1, 'complex_defects_var') and self.tab1.complex_defects_var.get():
            self.add_tab5()     ### Ajout pour restaurer l'onglet 5 si coché
        else:
            self.remove_tab5()  ### Ajout pour retirer l'onglet 5 si décoché
     
        ### Fin du bloc à ajouter

    def get_all_data(self):
        return {
            "tab1": self.tab1.get_data(),
            "tab2": self.tab2.get_data(),
            "tab3": self.tab3.get_data(),
        }

    def set_all_data(self, data):
        if "tab1" in data:
            self.tab1.set_data(data["tab1"])
        if "tab2" in data:
            self.tab2.set_data(data["tab2"])
        if "tab3" in data:
            self.tab3.set_data(data["tab3"])

    def save_file(self):
        data = self.get_all_data()
        file_path = filedialog.asksaveasfilename(defaultextension=".json",
                                                 filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if file_path:
            with open(file_path, "w") as f:
                json.dump(data, f, indent=2)
            messagebox.showinfo(self.lang_dict['save'], "Paramètres sauvegardés.")

    def open_file(self):
        file_path = filedialog.askopenfilename(defaultextension=".json",
                                               filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if file_path:
            with open(file_path, "r") as f:
                data = json.load(f)
            self.set_all_data(data)
            messagebox.showinfo(self.lang_dict['open'], "Paramètres chargés.")

    def open_fortran_input(self):
        file_path = filedialog.askopenfilename(
            defaultextension=".dat",
            filetypes=[("Input files", "*.dat;*.inp;*.txt"), ("All files", "*.*")]
        )
        if file_path:
            with open(file_path, "r") as f:
                lines = f.readlines()
            try:
                tab1_data, tab2_data, tab3_data = self.parse_fortran_input(lines)
                self.tab1.set_data(tab1_data)
                self.tab2.set_data(tab2_data)
                self.tab3.set_data(tab3_data)
                messagebox.showinfo("Import", "Paramètres importés depuis le fichier Fortran.")
            except Exception as e:
                messagebox.showerror("Erreur import", f"Erreur lors de la lecture : {e}")

    def parse_fortran_input(self, lines):
        """
        Parse un fichier d'entrée Fortran (format simplifié) et retourne un tuple (tab1_data, tab2_data, tab3_data)
        Cette version est minimaliste et doit être adaptée à la structure exacte du fichier!
        """
        tab1 = {}
        tab2 = {}
        tab3 = {}

        # Regex et helpers pour parser plus proprement si besoin
        import re

        # Etat du parseur
        section = None
        values = []
        raw = [l.strip() for l in lines if l.strip() and not l.strip().startswith("#") and not l.strip().startswith("-")]

        # --- Tab1 général ---
        # Ici, on suppose que les premiers paramètres sont dans l'ordre (exemple simpliste !)
        # Adapter selon la structure de tes fichiers réels
        try:
            # Commentaire = lignes jusqu'à la première valeur non textuelle
            comment_lines = []
            idx = 0
            while not re.match(r'^\d+(\.\d+)?$', raw[idx]):
                comment_lines.append(raw[idx])
                idx += 1
            tab1["comment"] = "\n".join(comment_lines)
            # Basename (nom du fichier de sortie)
            tab1["basename"] = raw[idx]
            idx += 1
            # Types chimiques
            tab1["chem_types_total"] = int(raw[idx])
            idx += 1
            tab1["chem_types_intrinsic"] = int(raw[idx])
            idx += 1
            # Sous-réseaux
            tab1["n_sublattices"] = int(raw[idx])
            idx += 1
            # ...à compléter selon l'ordre réel de tes paramètres...
        except Exception:
            pass  # Pour l'exemple, on ne plante pas si on n'a pas tout

        # --- Tab2 ---
        # Adapter ici selon le découpage réel, exemple :
        try:
            i = 0
            for idx, line in enumerate(raw):
                if "enrichissement" in line:
                    tab2["elem_enrichi"] = int(raw[idx+1])
                if "domaine de fraction" in line:
                    tab2["frac_min"] = float(raw[idx+1].split()[0])
                    tab2["frac_max"] = float(raw[idx+1].split()[1])
                if "référence pour les potentiels chimiques" in line:
                    vals = raw[idx+1].split()
                    tab2["type_ref"] = int(vals[0])
                    tab2["precision"] = float(vals[1].replace('D','E'))
                    tab2["niter"] = int(vals[2])
                if "Demi-largeurs des fenêtres" in line:
                    vals = raw[idx+1].split()
                    tab2["demi_largeurs"] = [str(v) for v in vals]
                if "Valeur initiale" in line:
                    vals = raw[idx+1].split()
                    tab2["mu2_init"] = float(vals[0])
                    tab2["mu2_nb"] = int(vals[1])
                    tab2["mu2_step"] = float(vals[2])
        except Exception:
            pass

        # --- Tab3 ---
        try:
            for idx, line in enumerate(raw):
                if "Energie de référence de la cellule sans défaut" in line:
                    tab3["energy_ref"] = float(raw[idx+1])
                if "Volume de référence" in line:
                    tab3["vol_ref"] = float(raw[idx+1])
                if "Nombre de mailles" in line:
                    tab3["nb_cells"] = int(raw[idx+1])
                if "Lacunes L(r)" in line and "brutes" in raw[idx-1]:
                    tab3["energy_lacune"] = float(raw[idx+1])
                if "Interstitiels 1(r)" in line and "brutes" in raw[idx-1]:
                    vals = raw[idx+1].split()
                    tab3["energy_inter1"] = [float(v) for v in vals]
                if "Substitutionnels 2(r)" in line and "brutes" in raw[idx-1]:
                    tab3["energy_subst2"] = float(raw[idx+1])
                if "Interstitiels 2(r)" in line and "brutes" in raw[idx-1]:
                    vals = raw[idx+1].split()
                    tab3["energy_inter2"] = [float(v) for v in vals]
                if "Lacunes L(r)" in line and "Volumes" in raw[idx-1]:
                    tab3["vol_lacune"] = float(raw[idx+1])
                if "Interstitiels 1(r)" in line and "Volumes" in raw[idx-1]:
                    vals = raw[idx+1].split()
                    tab3["vol_inter1"] = [float(v) for v in vals]
                if "Substitutionnels 2(r)" in line and "Volumes" in raw[idx-1]:
                    tab3["vol_subst2"] = float(raw[idx+1])
                if "Interstitiels 2(r)" in line and "Volumes" in raw[idx-1]:
                    vals = raw[idx+1].split()
                    tab3["vol_inter2"] = [float(v) for v in vals]
        except Exception:
            pass

        # Remplissage des valeurs par défaut si manquantes (sécurité)
        def_tab3 = Tab3Supercell(None, self.lang_dict).get_data()
        for k, v in def_tab3.items():
            if k not in tab3:
                tab3[k] = v
        def_tab2 = Tab2MuVTNPT(None, self.lang_dict, lambda: {"chem_types_total":2}).get_data()
        for k, v in def_tab2.items():
            if k not in tab2:
                tab2[k] = v
        def_tab1 = Tab1GeneralParameters(None, self.lang_dict).get_data()
        for k, v in def_tab1.items():
            if k not in tab1:
                tab1[k] = v

        return tab1, tab2, tab3

if __name__ == "__main__":
    app = FortranInputEditor()
    app.mainloop()
