import tkinter as tk
from tkinter import ttk

class Tab2MuVTNPT(ttk.Frame):
    """
    Onglet 2 : Paramètres muVT/NPT.
    Affiche dynamiquement les champs selon le mode de calcul et le nombre de types chimiques.
    """
    def __init__(self, parent, lang_dict, get_global_state):
        """
        get_global_state : fonction callback pour récupérer les infos globales (comme N_TYP, mode de calcul, etc)
        """
        super().__init__(parent)
        self.lang_dict = lang_dict
        self.get_global_state = get_global_state
        self.demi_largeur_entries = []
        self.create_widgets()

    def create_widgets(self):
        row = 0

        # Section 1 : Élément à enrichir
        ttk.Label(self, text="Élément à enrichir (index)").grid(row=row, column=0, sticky="w")
        self.elem_enrichi_var = tk.IntVar(value=2)
        self.elem_enrichi_entry = ttk.Entry(self, textvariable=self.elem_enrichi_var, width=5)
        self.elem_enrichi_entry.grid(row=row, column=1, sticky="w")
        row += 1

        # Section 2 : Domaine de fraction atomique
        ttk.Label(self, text="Fraction atomique min").grid(row=row, column=0, sticky="w")
        self.frac_min_var = tk.DoubleVar(value=0.0)
        ttk.Entry(self, textvariable=self.frac_min_var, width=8).grid(row=row, column=1, sticky="w")
        ttk.Label(self, text="Fraction atomique max").grid(row=row, column=2, sticky="w")
        self.frac_max_var = tk.DoubleVar(value=0.05)
        ttk.Entry(self, textvariable=self.frac_max_var, width=8).grid(row=row, column=3, sticky="w")
        row += 1

        # Section 3 : Spécifiques muVT
        ttk.Label(self, text="Type intrinsèque de référence").grid(row=row, column=0, sticky="w")
        self.type_ref_var = tk.IntVar(value=1)
        ttk.Entry(self, textvariable=self.type_ref_var, width=5).grid(row=row, column=1, sticky="w")
        ttk.Label(self, text="Précision boucle (POT_1)").grid(row=row, column=2, sticky="w")
        self.precision_var = tk.DoubleVar(value=1e-8)
        ttk.Entry(self, textvariable=self.precision_var, width=12).grid(row=row, column=3, sticky="w")
        ttk.Label(self, text="N° max itérations").grid(row=row, column=4, sticky="w")
        self.niter_var = tk.IntVar(value=1)
        ttk.Entry(self, textvariable=self.niter_var, width=5).grid(row=row, column=5, sticky="w")
        row += 1

        # Section 4 : Demi-largeurs (pour N_TYP > 2)
        ttk.Label(self, text="Demi-largeurs fenêtres (hors enrichi)").grid(row=row, column=0, sticky="w")
        self.demi_largeur_frame = ttk.Frame(self)
        self.demi_largeur_frame.grid(row=row, column=1, columnspan=5, sticky="w")
        self.update_demi_largeur_fields()
        row += 1

        # Section 5 : Série mu(addition 2)
        ttk.Separator(self, orient='horizontal').grid(row=row, column=0, columnspan=6, sticky="we", pady=8)
        row += 1
        ttk.Label(self, text="Série mu(addition 2) :").grid(row=row, column=0, sticky="w")
        row += 1
        ttk.Label(self, text="Valeur initiale (eV)").grid(row=row, column=0, sticky="w")
        self.mu2_init_var = tk.DoubleVar(value=-12.0)
        ttk.Entry(self, textvariable=self.mu2_init_var, width=8).grid(row=row, column=1, sticky="w")
        ttk.Label(self, text="Nombre de valeurs").grid(row=row, column=2, sticky="w")
        self.mu2_nb_var = tk.IntVar(value=300)
        ttk.Entry(self, textvariable=self.mu2_nb_var, width=6).grid(row=row, column=3, sticky="w")
        ttk.Label(self, text="Pas (eV)").grid(row=row, column=4, sticky="w")
        self.mu2_step_var = tk.DoubleVar(value=0.01)
        ttk.Entry(self, textvariable=self.mu2_step_var, width=8).grid(row=row, column=5, sticky="w")

        # Pour la mise à jour dynamique si le nombre de types chimiques change
        self.after(500, self.check_update_fields)

    def get_n_types(self):
        # Récupère le nombre d'espèces chimiques depuis l'extérieur
        state = self.get_global_state()
        return state.get("chem_types_total", 2)

    def update_demi_largeur_fields(self):
        # Met à jour le nombre de champs de demi-largeur selon N_TYP
        n_types = self.get_n_types()
        n_fields = max(n_types - 1, 0)
        old_values = [v.get() for v in self.demi_largeur_entries]
        # Clear old widgets
        for widget in self.demi_largeur_frame.winfo_children():
            widget.destroy()
        self.demi_largeur_entries = []
        for i in range(n_fields):
            var = tk.StringVar()
            if i < len(old_values):
                var.set(old_values[i])
            else:
                var.set("0.00005")
            entry = ttk.Entry(self.demi_largeur_frame, textvariable=var, width=10)
            entry.grid(row=0, column=i)
            self.demi_largeur_entries.append(var)

    def check_update_fields(self):
        # Vérifie si le nombre de types chimiques a changé, et met à jour si besoin
        self.update_demi_largeur_fields()
        self.after(500, self.check_update_fields)

    def get_data(self):
        return {
            "elem_enrichi": self.elem_enrichi_var.get(),
            "frac_min": self.frac_min_var.get(),
            "frac_max": self.frac_max_var.get(),
            "type_ref": self.type_ref_var.get(),
            "precision": self.precision_var.get(),
            "niter": self.niter_var.get(),
            "demi_largeurs": [v.get() for v in self.demi_largeur_entries],
            "mu2_init": self.mu2_init_var.get(),
            "mu2_nb": self.mu2_nb_var.get(),
            "mu2_step": self.mu2_step_var.get(),
        }

    def set_data(self, data):
        self.elem_enrichi_var.set(data.get("elem_enrichi", 2))
        self.frac_min_var.set(data.get("frac_min", 0.0))
        self.frac_max_var.set(data.get("frac_max", 0.05))
        self.type_ref_var.set(data.get("type_ref", 1))
        self.precision_var.set(data.get("precision", 1e-8))
        self.niter_var.set(data.get("niter", 1))
        # Champs dynamiques
        self.update_demi_largeur_fields()
        demi = data.get("demi_largeurs", [])
        for i, var in enumerate(self.demi_largeur_entries):
            if i < len(demi):
                var.set(demi[i])
        self.mu2_init_var.set(data.get("mu2_init", -12.0))
        self.mu2_nb_var.set(data.get("mu2_nb", 300))
        self.mu2_step_var.set(data.get("mu2_step", 0.01))

    def set_labels(self, lang_dict):
        self.lang_dict = lang_dict
        for widget in self.winfo_children():
            widget.destroy()
        self.create_widgets()
