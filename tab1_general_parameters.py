import tkinter as tk
from tkinter import ttk

class Tab1GeneralParameters(ttk.Frame):
    """
    Onglet 1 : Paramètres généraux du fichier d'entrée.
    Champs dynamiques pour 'sites par maille', 'énergies de référence', 'sous-réseaux occupés', etc.
    """
    def __init__(self, parent, lang_dict, tab_manager=None):
        super().__init__(parent)
        self.lang_dict = lang_dict
        self.tab_manager = tab_manager
        self.sous_reseaux_occ_vars = []    ### Ajout pour le champ dynamique "sous-réseaux occupés"
        self.sites_per_cell_entries = []
        self.ref_energies_entries = []
        self.create_widgets()

    def create_widgets(self):
        row = 0
        # Commentaire (multiligne)
        default_comment = (
            "Ti hcp avec O en  substitution et  en insertion\n"
            "(sites inter. = octa + tetra specifies dans cet ordre)\n"
            "Ti = 1 - O = 2"
        )
        ttk.Label(self, text=self.lang_dict['comment']).grid(row=row, column=0, sticky="w")
        self.comment_entry = tk.Text(self, height=3, width=60)
        self.comment_entry.grid(row=row, column=1, columnspan=4, sticky="we")
        self.comment_entry.insert("1.0", default_comment)
        row += 1

        # Base du nom des fichiers de sortie
        ttk.Label(self, text=self.lang_dict['basename']).grid(row=row, column=0, sticky="w")
        self.basename_var = tk.StringVar(value="Ti_1000K")
        ttk.Entry(self, textvariable=self.basename_var, width=30).grid(row=row, column=1, columnspan=4, sticky="w")
        row += 1

        # Nombres de types chimiques (total, intrinsèque)
        ttk.Label(self, text=self.lang_dict['chem_types']).grid(row=row, column=0, sticky="w")
        self.chem_types_total_var = tk.IntVar(value=2)
        self.chem_types_intrinsic_var = tk.IntVar(value=1)
        chem_total_entry = ttk.Entry(self, textvariable=self.chem_types_total_var, width=5)
        chem_intrin_entry = ttk.Entry(self, textvariable=self.chem_types_intrinsic_var, width=5)
        chem_total_entry.grid(row=row, column=1, sticky="w")
        chem_intrin_entry.grid(row=row, column=2, sticky="w")
        chem_total_entry.bind("<FocusOut>", lambda e: self.update_ref_energies_fields())
        chem_total_entry.bind("<Return>", lambda e: self.update_ref_energies_fields())
        chem_intrin_entry.bind("<FocusOut>", lambda e: [self.update_ref_energies_fields(), self.update_sous_reseaux_occ_fields()])   ### Ajout update_sous_reseaux_occ_fields
        chem_intrin_entry.bind("<Return>", lambda e: [self.update_ref_energies_fields(), self.update_sous_reseaux_occ_fields()])     ### Ajout update_sous_reseaux_occ_fields
        row += 1

        # Nombre de sous-réseaux
        ttk.Label(self, text=self.lang_dict['n_sublattices']).grid(row=row, column=0, sticky="w")
        self.n_sublattices_var = tk.IntVar(value=3)
        n_sublattices_entry = ttk.Entry(self, textvariable=self.n_sublattices_var, width=5)
        n_sublattices_entry.grid(row=row, column=1, sticky="w")
        n_sublattices_entry.bind("<FocusOut>", lambda e: self.update_sites_per_cell_fields())
        n_sublattices_entry.bind("<Return>", lambda e: self.update_sites_per_cell_fields())
        row += 1

        # Présence de sous-réseaux interstitiels (O/o ou N/n)
        ttk.Label(self, text=self.lang_dict['interstitial']).grid(row=row, column=0, sticky="w")
        self.interstitial_var = tk.BooleanVar()
        ttk.Checkbutton(self, variable=self.interstitial_var, text=self.lang_dict['interstitial_present']).grid(row=row, column=1, sticky="w")
        row += 1

        # Prise en compte des interstitiels intrinsèques (O/o ou N/n)
        ttk.Label(self, text=self.lang_dict['intrinsic_interstitial']).grid(row=row, column=0, sticky="w")
        self.intrinsic_interstitial_var = tk.BooleanVar()
        ttk.Checkbutton(self, variable=self.intrinsic_interstitial_var, text=self.lang_dict['intrinsic_interstitial_present']).grid(row=row, column=1, sticky="w")
        row += 1

        # Nombres de sous-réseaux occupés (état fondamental) - DYNAMIQUE
        ttk.Label(self, text=self.lang_dict['n_sublattices_intrinsic']).grid(row=row, column=0, sticky="w")
        self.sous_reseaux_occ_frame = ttk.Frame(self)    ### Ajout d'un frame dynamique
        self.sous_reseaux_occ_frame.grid(row=row, column=1, columnspan=4, sticky="w")
        self.update_sous_reseaux_occ_fields() ### Ajout appel initial
        row += 1

        # Nombre de sites par maille pour chaque sous-réseau (dynamique)
        ttk.Label(self, text=self.lang_dict['sites_per_cell']).grid(row=row, column=0, sticky="w")
        self.sites_per_cell_frame = ttk.Frame(self)
        self.sites_per_cell_frame.grid(row=row, column=1, columnspan=4, sticky="w")
        self.update_sites_per_cell_fields(default_values=["2", "2", "4"])
        row += 1

        # Défauts complexes (O/o ou N/n)
        ttk.Label(self, text=self.lang_dict['complex_defects']).grid(row=row, column=0, sticky="w")
        self.complex_defects_var = tk.BooleanVar()
        ttk.Checkbutton(self, variable=self.complex_defects_var, text=self.lang_dict['complex_defects_present']).grid(row=row, column=1, sticky="w")
        self.complex_defects_var.trace_add("write", lambda *args: self.on_complex_defects_changed())
        row += 1

        # Température (K) par défaut à 1000.0
        ttk.Label(self, text=self.lang_dict['temperature']).grid(row=row, column=0, sticky="w")
        self.temperature_var = tk.DoubleVar(value=1000.0)
        ttk.Entry(self, textvariable=self.temperature_var, width=8).grid(row=row, column=1, sticky="w")
        row += 1

        # Pression externe (kbar)
        ttk.Label(self, text=self.lang_dict['pressure']).grid(row=row, column=0, sticky="w")
        self.pressure_var = tk.DoubleVar()
        ttk.Entry(self, textvariable=self.pressure_var, width=8).grid(row=row, column=1, sticky="w")
        row += 1

        # Énergies de référence (dynamique, côte à côte)
        ttk.Label(self, text=self.lang_dict['ref_energies']).grid(row=row, column=0, sticky="w")
        self.ref_energies_frame = ttk.Frame(self)
        self.ref_energies_frame.grid(row=row, column=1, columnspan=4, sticky="w")
        self.update_ref_energies_fields(default_values=["0.00", "0.00"])
        row += 1

        # Grandeurs thermodynamiques (A/a ou M/m) côte à côte
        ttk.Label(self, text=self.lang_dict['thermo_mode']).grid(row=row, column=0, sticky="w")
        self.thermo_mode_var = tk.StringVar(value="a")
        ttk.Radiobutton(self, text=self.lang_dict['by_atom'], variable=self.thermo_mode_var, value="a").grid(row=row, column=1, sticky="w")
        ttk.Radiobutton(self, text=self.lang_dict['by_cell'], variable=self.thermo_mode_var, value="m").grid(row=row, column=2, sticky="w")
        row += 1

        # Indicateur d'écriture de l'énergie libre par atome (T/t ou F/f) côte à côte
        ttk.Label(self, text=self.lang_dict['free_energy_type']).grid(row=row, column=0, sticky="w")
        self.free_energy_type_var = tk.StringVar(value="t")
        ttk.Radiobutton(self, text=self.lang_dict['total'], variable=self.free_energy_type_var, value="t").grid(row=row, column=1, sticky="w")
        ttk.Radiobutton(self, text=self.lang_dict['formation'], variable=self.free_energy_type_var, value="f").grid(row=row, column=2, sticky="w")
        row += 1

        # Type de calcul (muVT, NPT, NPT=0)
        ttk.Label(self, text=self.lang_dict['calc_type']).grid(row=row, column=0, sticky="w")
        self.calc_type_var = tk.StringVar(value="muVT")
        ttk.Combobox(self, textvariable=self.calc_type_var, values=["muVT", "NPT", "NPT=0"], width=8).grid(row=row, column=1, sticky="w")
        row += 1

        # Défauts ponctuels chargés (Checkbutton 0/1)
        ttk.Label(self, text="Défauts ponctuels chargés").grid(row=row, column=0, sticky="w")   ### Label corrigé
        self.defaut_charge_var = tk.IntVar(value=0)      ### Ajout
        ttk.Checkbutton(self, variable=self.defaut_charge_var, command=self.on_defaut_charge_changed).grid(row=row, column=1, sticky="w")   ### Ajout
        row += 1

        self.grid_columnconfigure(1, weight=1)

    def on_complex_defects_changed(self):
        if self.tab_manager:
            if self.complex_defects_var.get():
                self.tab_manager.add_tab5()
            else:
                self.tab_manager.remove_tab5()

    def update_sous_reseaux_occ_fields(self):
        # Crée autant de cases que le nombre de types chimiques intrinsèques
        for widget in self.sous_reseaux_occ_frame.winfo_children():
            widget.destroy()
        self.sous_reseaux_occ_vars = []
        try:
            n = int(self.chem_types_intrinsic_var.get())
        except Exception:
            n = 0
        for i in range(n):
            var = tk.IntVar(value=1)
            entry = ttk.Entry(self.sous_reseaux_occ_frame, textvariable=var, width=5)
            entry.grid(row=0, column=i)
            self.sous_reseaux_occ_vars.append(var)

    def update_sites_per_cell_fields(self, default_values=None):
        old_values = [entry.get() for entry in self.sites_per_cell_entries]
        for widget in self.sites_per_cell_frame.winfo_children():
            widget.destroy()
        self.sites_per_cell_entries = []
        try:
            n = int(self.n_sublattices_var.get())
        except Exception:
            n = 0
        for i in range(n):
            var = tk.StringVar()
            if default_values and i < len(default_values):
                var.set(default_values[i])
            elif i < len(old_values):
                var.set(old_values[i])
            entry = ttk.Entry(self.sites_per_cell_frame, textvariable=var, width=5)
            entry.grid(row=0, column=i)
            self.sites_per_cell_entries.append(var)

    def update_ref_energies_fields(self, default_values=None):
        old_values = [entry.get() for entry in self.ref_energies_entries]
        for widget in self.ref_energies_frame.winfo_children():
            widget.destroy()
        self.ref_energies_entries = []
        try:
            n = int(self.chem_types_total_var.get())
        except Exception:
            n = 0
        for i in range(n):
            var = tk.StringVar()
            if default_values and i < len(default_values):
                var.set(default_values[i])
            elif i < len(old_values):
                var.set(old_values[i])
            entry = ttk.Entry(self.ref_energies_frame, textvariable=var, width=7)
            entry.grid(row=0, column=i)
            self.ref_energies_entries.append(var)

    def on_defaut_charge_changed(self):
        if self.tab_manager:
            if self.defaut_charge_var.get():
                self.tab_manager.add_tab4()
            else:
                self.tab_manager.remove_tab4()

    def get_data(self):
        return {
            "comment": self.comment_entry.get("1.0", tk.END).strip(),
            "basename": self.basename_var.get(),
            "chem_types_total": self.chem_types_total_var.get(),
            "chem_types_intrinsic": self.chem_types_intrinsic_var.get(),
            "n_sublattices": self.n_sublattices_var.get(),
            "interstitial": self.interstitial_var.get(),
            "intrinsic_interstitial": self.intrinsic_interstitial_var.get(),
            "n_sublattices_intrinsic": [v.get() for v in self.sous_reseaux_occ_vars],   ### Nouvelle clé : liste dynamique
            "sites_per_cell": [v.get() for v in self.sites_per_cell_entries],
            "complex_defects": self.complex_defects_var.get(),
            "temperature": self.temperature_var.get(),
            "pressure": self.pressure_var.get(),
            "ref_energies": [v.get() for v in self.ref_energies_entries],
            "thermo_mode": self.thermo_mode_var.get(),
            "free_energy_type": self.free_energy_type_var.get(),
            "calc_type": self.calc_type_var.get(),
            "defaut_charge": self.defaut_charge_var.get(),    ### Clé corrigée (au lieu de dp_charges)
        }

    def set_data(self, data):
        self.comment_entry.delete("1.0", tk.END)
        self.comment_entry.insert("1.0", data.get("comment", ""))
        self.basename_var.set(data.get("basename", "Ti_1000K"))
        self.chem_types_total_var.set(data.get("chem_types_total", 2))
        self.chem_types_intrinsic_var.set(data.get("chem_types_intrinsic", 1))
        self.n_sublattices_var.set(data.get("n_sublattices", 3))
        self.interstitial_var.set(data.get("interstitial", False))
        self.intrinsic_interstitial_var.set(data.get("intrinsic_interstitial", False))
        self.update_sous_reseaux_occ_fields()
        sous_reseaux = data.get("n_sublattices_intrinsic", [1]*self.chem_types_intrinsic_var.get())
        for i, var in enumerate(self.sous_reseaux_occ_vars):
            if i < len(sous_reseaux):
                var.set(sous_reseaux[i])
        self.update_sites_per_cell_fields()
        sites = data.get("sites_per_cell", [])
        for i, var in enumerate(self.sites_per_cell_entries):
            if i < len(sites):
                var.set(sites[i])
        self.complex_defects_var.set(data.get("complex_defects", False))
        self.temperature_var.set(data.get("temperature", 1000.0))
        self.pressure_var.set(data.get("pressure", 0.0))
        self.update_ref_energies_fields()
        energies = data.get("ref_energies", [])
        for i, var in enumerate(self.ref_energies_entries):
            if i < len(energies):
                var.set(energies[i])
        self.thermo_mode_var.set(data.get("thermo_mode", "a"))
        self.free_energy_type_var.set(data.get("free_energy_type", "t"))
        self.calc_type_var.set(data.get("calc_type", "muVT"))
        self.defaut_charge_var.set(data.get("defaut_charge", 0))    ### Correction

    def set_labels(self, lang_dict):
        self.lang_dict = lang_dict
        for widget in self.winfo_children():
            widget.destroy()
        self.create_widgets()
