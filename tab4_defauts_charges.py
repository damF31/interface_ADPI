import tkinter as tk
from tkinter import ttk

class Tab4DefautsCharges(ttk.Frame):
    """Onglet 4 : Options défauts ponctuels chargés"""

    def __init__(self, parent, lang_dict=None):
        super().__init__(parent)
        self.lang_dict = lang_dict or {}
        self.create_widgets()

    def create_widgets(self):
        row = 0
        label = self.lang_dict.get('tab4', "Options défauts chargés")
        ttk.Label(self, text=label, font=("Arial", 12, "bold")).grid(row=row, column=0, sticky="w", pady=10)
        row += 1

        # Exemples de champs spécifiques à cet onglet (à adapter à tes besoins)
        ttk.Label(self, text="Exemple : Niveau de Fermi (eV)").grid(row=row, column=0, sticky="w")
        self.fermi_level_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.fermi_level_var, width=10).grid(row=row, column=1, sticky="w")
        row += 1

        ttk.Label(self, text="Exemple : Correction image charge").grid(row=row, column=0, sticky="w")
        self.correction_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(self, variable=self.correction_var, text="Correction activée").grid(row=row, column=1, sticky="w")
        row += 1

    def get_data(self):
        return {
            "fermi_level": self.fermi_level_var.get(),
            "correction_active": self.correction_var.get()
        }

    def set_data(self, data):
        self.fermi_level_var.set(data.get("fermi_level", 0.0))
        self.correction_var.set(data.get("correction_active", False))
