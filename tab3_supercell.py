import tkinter as tk
from tkinter import ttk

class Tab3Supercell(ttk.Frame):
    """
    Onglet 3 : Paramètres de la supercellule de référence et énergies/volumes des défauts ponctuels.
    """
    def __init__(self, parent, lang_dict):
        super().__init__(parent)
        self.lang_dict = lang_dict
        self.create_widgets()

    def create_widgets(self):
        row = 0
        # --- Paramètres cellule de référence ---
        title1 = ttk.Label(self, text="Paramètres de la supercellule de référence", font=("Arial", 10, "bold"))
        title1.grid(row=row, column=0, columnspan=4, sticky="w", pady=(10, 2))
        row += 1

        ttk.Label(self, text="Énergie de référence (eV)").grid(row=row, column=0, sticky="w")
        self.energy_ref_var = tk.DoubleVar(value=-749.07575481)
        ttk.Entry(self, textvariable=self.energy_ref_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        ttk.Label(self, text="Volume de référence (Å³)").grid(row=row, column=0, sticky="w")
        self.vol_ref_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.vol_ref_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        ttk.Label(self, text="Nombre de mailles").grid(row=row, column=0, sticky="w")
        self.nb_cells_var = tk.IntVar(value=48)
        ttk.Entry(self, textvariable=self.nb_cells_var, width=8).grid(row=row, column=1, sticky="w")
        row += 2

        # --- Energies brutes des défauts ponctuels ---
        title2 = ttk.Label(self, text="Énergies 'brutes' des défauts ponctuels (eV)", font=("Arial", 10, "bold"))
        title2.grid(row=row, column=0, columnspan=4, sticky="w", pady=(10, 2))
        row += 1

        # Lacunes L(r)
        ttk.Label(self, text="Lacunes L(r)").grid(row=row, column=0, sticky="w")
        self.energy_lacune_var = tk.DoubleVar(value=-739.2145188)
        ttk.Entry(self, textvariable=self.energy_lacune_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        # Interstitiels 1(r)
        ttk.Label(self, text="Interstitiels 1(r)").grid(row=row, column=0, sticky="w")
        self.energy_inter1_var1 = tk.DoubleVar(value=-700.236)
        self.energy_inter1_var2 = tk.DoubleVar(value=-700.236)
        ttk.Entry(self, textvariable=self.energy_inter1_var1, width=16).grid(row=row, column=1, sticky="w")
        ttk.Entry(self, textvariable=self.energy_inter1_var2, width=16).grid(row=row, column=2, sticky="w")
        row += 1

        # Substitutionnels 2(r)
        ttk.Label(self, text="Substitutionnels 2(r)").grid(row=row, column=0, sticky="w")
        self.energy_subst2_var = tk.DoubleVar(value=-749.6382079)
        ttk.Entry(self, textvariable=self.energy_subst2_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        # Interstitiels 2(r)
        ttk.Label(self, text="Interstitiels 2(r)").grid(row=row, column=0, sticky="w")
        self.energy_inter2_var1 = tk.DoubleVar(value=-759.5010916)
        self.energy_inter2_var2 = tk.DoubleVar(value=-758.2618091)
        ttk.Entry(self, textvariable=self.energy_inter2_var1, width=16).grid(row=row, column=1, sticky="w")
        ttk.Entry(self, textvariable=self.energy_inter2_var2, width=16).grid(row=row, column=2, sticky="w")
        row += 2

        # --- Volumes bruts des défauts ponctuels ---
        title3 = ttk.Label(self, text="Volumes 'bruts' des défauts ponctuels (Å³)", font=("Arial", 10, "bold"))
        title3.grid(row=row, column=0, columnspan=4, sticky="w", pady=(10, 2))
        row += 1

        # Lacunes L(r)
        ttk.Label(self, text="Lacunes L(r)").grid(row=row, column=0, sticky="w")
        self.vol_lacune_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.vol_lacune_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        # Interstitiels 1(r)
        ttk.Label(self, text="Interstitiels 1(r)").grid(row=row, column=0, sticky="w")
        self.vol_inter1_var1 = tk.DoubleVar(value=0.0)
        self.vol_inter1_var2 = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.vol_inter1_var1, width=16).grid(row=row, column=1, sticky="w")
        ttk.Entry(self, textvariable=self.vol_inter1_var2, width=16).grid(row=row, column=2, sticky="w")
        row += 1

        # Substitutionnels 2(r)
        ttk.Label(self, text="Substitutionnels 2(r)").grid(row=row, column=0, sticky="w")
        self.vol_subst2_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.vol_subst2_var, width=16).grid(row=row, column=1, sticky="w")
        row += 1

        # Interstitiels 2(r)
        ttk.Label(self, text="Interstitiels 2(r)").grid(row=row, column=0, sticky="w")
        self.vol_inter2_var1 = tk.DoubleVar(value=0.0)
        self.vol_inter2_var2 = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.vol_inter2_var1, width=16).grid(row=row, column=1, sticky="w")
        ttk.Entry(self, textvariable=self.vol_inter2_var2, width=16).grid(row=row, column=2, sticky="w")
        row += 1

        self.grid_columnconfigure(1, weight=1)

    def get_data(self):
        return {
            "energy_ref": self.energy_ref_var.get(),
            "vol_ref": self.vol_ref_var.get(),
            "nb_cells": self.nb_cells_var.get(),
            "energy_lacune": self.energy_lacune_var.get(),
            "energy_inter1": [self.energy_inter1_var1.get(), self.energy_inter1_var2.get()],
            "energy_subst2": self.energy_subst2_var.get(),
            "energy_inter2": [self.energy_inter2_var1.get(), self.energy_inter2_var2.get()],
            "vol_lacune": self.vol_lacune_var.get(),
            "vol_inter1": [self.vol_inter1_var1.get(), self.vol_inter1_var2.get()],
            "vol_subst2": self.vol_subst2_var.get(),
            "vol_inter2": [self.vol_inter2_var1.get(), self.vol_inter2_var2.get()],
        }

    def set_data(self, data):
        self.energy_ref_var.set(data.get("energy_ref", -749.07575481))
        self.vol_ref_var.set(data.get("vol_ref", 0.0))
        self.nb_cells_var.set(data.get("nb_cells", 48))
        self.energy_lacune_var.set(data.get("energy_lacune", -739.2145188))
        energy_inter1 = data.get("energy_inter1", [-700.236, -700.236])
        self.energy_inter1_var1.set(energy_inter1[0])
        self.energy_inter1_var2.set(energy_inter1[1])
        self.energy_subst2_var.set(data.get("energy_subst2", -749.6382079))
        energy_inter2 = data.get("energy_inter2", [-759.5010916, -758.2618091])
        self.energy_inter2_var1.set(energy_inter2[0])
        self.energy_inter2_var2.set(energy_inter2[1])
        self.vol_lacune_var.set(data.get("vol_lacune", 0.0))
        vol_inter1 = data.get("vol_inter1", [0.0, 0.0])
        self.vol_inter1_var1.set(vol_inter1[0])
        self.vol_inter1_var2.set(vol_inter1[1])
        self.vol_subst2_var.set(data.get("vol_subst2", 0.0))
        vol_inter2 = data.get("vol_inter2", [0.0, 0.0])
        self.vol_inter2_var1.set(vol_inter2[0])
        self.vol_inter2_var2.set(vol_inter2[1])

    def set_labels(self, lang_dict):
        self.lang_dict = lang_dict
        for widget in self.winfo_children():
            widget.destroy()
        self.create_widgets()
