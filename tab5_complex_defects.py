import tkinter as tk
from tkinter import ttk

class ComplexPanel(ttk.LabelFrame):
    """Affiche un complexe avec ses paramètres dans l'onglet 5."""
    def __init__(self, parent, index):
        super().__init__(parent, text=f"Complexe #{index+1}")
        self.index = index
        self.create_widgets()

    def create_widgets(self):
        row = 0
        self.num_var = tk.IntVar(value=self.index + 1)
        self.multip_var = tk.IntVar(value=1)
        self.sr_corresp_var = tk.IntVar(value=1)
        self.n_sites_var = tk.IntVar(value=2)
        ttk.Label(self, text="Multiplicité d’orientation").grid(row=row, column=0, sticky="w")
        ttk.Entry(self, textvariable=self.multip_var, width=5).grid(row=row, column=1)
        ttk.Label(self, text="Sous-réseau correspondant").grid(row=row, column=2, sticky="w")
        ttk.Entry(self, textvariable=self.sr_corresp_var, width=5).grid(row=row, column=3)
        ttk.Label(self, text="Nombre de sites").grid(row=row, column=4, sticky="w")
        ttk.Entry(self, textvariable=self.n_sites_var, width=5).grid(row=row, column=5)
        row += 1

        # Champs dynamiques pour sous-réseaux occupés et types chimiques
        ttk.Label(self, text="Numéros de sous-réseaux occupés").grid(row=row, column=0, sticky="w")
        self.sr_sites_vars = []
        for i in range(5):  # max 5, ajustable
            var = tk.IntVar(value=1)
            entry = ttk.Entry(self, textvariable=var, width=4)
            entry.grid(row=row, column=i+1)
            self.sr_sites_vars.append(var)
        row += 1
        ttk.Label(self, text="Types chimiques sur ces sites").grid(row=row, column=0, sticky="w")
        self.types_sites_vars = []
        for i in range(5):
            var = tk.IntVar(value=0)
            entry = ttk.Entry(self, textvariable=var, width=4)
            entry.grid(row=row, column=i+1)
            self.types_sites_vars.append(var)
        row += 1

        # Energie et volume
        ttk.Label(self, text="Energie (eV)").grid(row=row, column=0, sticky="w")
        self.energy_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.energy_var, width=10).grid(row=row, column=1)
        ttk.Label(self, text="Volume (Å^3)").grid(row=row, column=2, sticky="w")
        self.volume_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.volume_var, width=10).grid(row=row, column=3)
        row += 1

    def get_data(self):
        return {
            "num": self.num_var.get(),
            "multiplicity": self.multip_var.get(),
            "sr_corresp": self.sr_corresp_var.get(),
            "n_sites": self.n_sites_var.get(),
            "sr_sites": [v.get() for v in self.sr_sites_vars],
            "types_sites": [v.get() for v in self.types_sites_vars],
            "energy": self.energy_var.get(),
            "volume": self.volume_var.get()
        }

    def set_data(self, data):
        self.num_var.set(data.get("num", self.index+1))
        self.multip_var.set(data.get("multiplicity", 1))
        self.sr_corresp_var.set(data.get("sr_corresp", 1))
        self.n_sites_var.set(data.get("n_sites", 2))
        sr_sites = data.get("sr_sites", [1]*len(self.sr_sites_vars))
        for i, var in enumerate(self.sr_sites_vars):
            if i < len(sr_sites):
                var.set(sr_sites[i])
        types_sites = data.get("types_sites", [0]*len(self.types_sites_vars))
        for i, var in enumerate(self.types_sites_vars):
            if i < len(types_sites):
                var.set(types_sites[i])
        self.energy_var.set(data.get("energy", 0.0))
        self.volume_var.set(data.get("volume", 0.0))


class Tab5ComplexDefects(ttk.Frame):
    def __init__(self, parent, lang_dict):
        super().__init__(parent)
        self.lang_dict = lang_dict
        self.panels = []
        self.create_widgets()

    def create_widgets(self):
        row = 0
        ttk.Label(self, text="Nombre maximum de sites des complexes").grid(row=row, column=0, sticky="w")
        self.max_sites_var = tk.IntVar(value=20)
        ttk.Entry(self, textvariable=self.max_sites_var, width=8).grid(row=row, column=1)
        row += 1

        ttk.Label(self, text="Nombre de types de défauts complexes").grid(row=row, column=0, sticky="w")
        self.n_complex_types_var = tk.IntVar(value=1)
        entry = ttk.Entry(self, textvariable=self.n_complex_types_var, width=5)
        entry.grid(row=row, column=1)
        entry.bind("<FocusOut>", lambda e: self.update_panels())
        entry.bind("<Return>", lambda e: self.update_panels())
        row += 1

        self.panels_frame = ttk.Frame(self)
        self.panels_frame.grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1

        self.update_panels()

    def update_panels(self):
        for widget in self.panels_frame.winfo_children():
            widget.destroy()
        self.panels = []
        n = self.n_complex_types_var.get()
        for i in range(n):
            panel = ComplexPanel(self.panels_frame, i)
            panel.grid(row=i, column=0, pady=10, sticky="w")
            self.panels.append(panel)

    def get_data(self):
        return {
            "max_sites": self.max_sites_var.get(),
            "n_complex_types": self.n_complex_types_var.get(),
            "complexes": [panel.get_data() for panel in self.panels]
        }

    def set_data(self, data):
        self.max_sites_var.set(data.get("max_sites", 20))
        self.n_complex_types_var.set(data.get("n_complex_types", 1))
        self.update_panels()
        complexes = data.get("complexes", [])
        for i, panel in enumerate(self.panels):
            if i < len(complexes):
                panel.set_data(complexes[i])
