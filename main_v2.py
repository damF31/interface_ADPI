import tkinter as tk
from tkinter import ttk, filedialog, messagebox
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
        filemenu.add_command(label="Exporter .dat/.inp...", command=self.export_fortran_input)
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
        # Correction: vérifier si self.tab5 existe et est bien un onglet du notebook
        if hasattr(self, 'tab5_exists') and self.tab5_exists and hasattr(self, 'tab5'):
            try:
                idx = self.tabs.index(self.tab5)
                self.tabs.forget(idx)
            except Exception:
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
        # Bloc de rétablissement des onglets dynamiques
        if hasattr(self.tab1, 'defaut_charge_var') and self.tab1.defaut_charge_var.get():
            self.add_tab4()
        else:
            self.remove_tab4()
        if hasattr(self.tab1, 'complex_defects_var') and self.tab1.complex_defects_var.get():
            self.add_tab5()
        else:
            self.remove_tab5()

    def get_all_data(self):
        data = {
            "tab1": self.tab1.get_data(),
            "tab2": self.tab2.get_data(),
            "tab3": self.tab3.get_data(),
        }
        # Ajoute ici tab4/tab5 si existants
        if hasattr(self, 'tab4_exists') and self.tab4_exists and hasattr(self, 'tab4'):
            if hasattr(self.tab4, 'get_data'):
                data['tab4'] = self.tab4.get_data()
        if hasattr(self, 'tab5_exists') and self.tab5_exists and hasattr(self, 'tab5'):
            if hasattr(self.tab5, 'get_data'):
                data['tab5'] = self.tab5.get_data()
        return data

    def set_all_data(self, data):
        if "tab1" in data:
            self.tab1.set_data(data["tab1"])
        if "tab2" in data:
            self.tab2.set_data(data["tab2"])
        if "tab3" in data:
            self.tab3.set_data(data["tab3"])
        if "tab4" in data and hasattr(self, 'tab4'):
            self.tab4.set_data(data["tab4"])
        if "tab5" in data and hasattr(self, 'tab5'):
            self.tab5.set_data(data["tab5"])

    # UTILISATION DU MODULE CENTRALISE file_io.py POUR LES OPERATIONS DE FICHIER

    def save_file(self):
        """Sauvegarde de la session au format JSON."""
        data = self.get_all_data()
        save_session(data, parent_window=self)

    def open_file(self):
        """Chargement d'une session au format JSON."""
        data = load_session(parent_window=self)
        if data:
            self.set_all_data(data)
            messagebox.showinfo(self.lang_dict['open'], "Paramètres chargés.")

    def export_fortran_input(self):
        """Export au format texte FORTRAN (.dat/.inp)."""
        data = self.get_all_data()
        export_fortran(data, parent_window=self)

    def open_fortran_input(self):
        """Import d’un fichier texte FORTRAN (.dat/.inp)."""
        data = import_fortran(parent_window=self)
        if data:
            # TODO: Selon le mapping, adapter la façon dont on répartit les données dans les onglets
            # Pour l’instant on suppose que les clés correspondent à tab1/tab2/tab3/tab4/tab5
            if "tab1" in data:
                self.tab1.set_data(data["tab1"])
            if "tab2" in data:
                self.tab2.set_data(data["tab2"])
            if "tab3" in data:
                self.tab3.set_data(data["tab3"])
            if hasattr(self, 'tab4') and "tab4" in data:
                self.tab4.set_data(data["tab4"])
            if hasattr(self, 'tab5') and "tab5" in data:
                self.tab5.set_data(data["tab5"])
            messagebox.showinfo("Import", "Paramètres importés depuis le fichier Fortran.")

if __name__ == "__main__":
    app = FortranInputEditor()
    app.mainloop()
