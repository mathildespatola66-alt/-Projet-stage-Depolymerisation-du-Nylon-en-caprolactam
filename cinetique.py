import os
import json
import numpy as np
import pandas as pd
import itertools
from scipy.optimize import curve_fit

# IMPORTATION DES OUTILS DE TON JUMEAU NUMÉRIQUE
from random_forest import charger_et_preparer_donnees, entrainer_expert, DOSSIER_TRAVAIL

# ==========================================
# 1. LE MODÈLE MATHÉMATIQUE RIGOUREUX
# ==========================================
def modele_cinetique_n_m(X_data, A, Ea, n, m):
    """ Équation différentielle intégrée : dX/dt = A*exp(-Ea/RT) * [cat]^m * (1-X)^n """
    t, T, cat = X_data
    R = 8.314
    
    # Éviter l'erreur 0^m
    cat_eff = np.where(cat <= 0, 1e-5, cat)
    k = A * np.exp(-Ea / (R * T)) * (cat_eff**m)
    
    X_pred = np.zeros_like(t)
    for i in range(len(t)):
        if abs(n - 1.0) < 1e-4:  # Ordre 1
            X_pred[i] = 1.0 - np.exp(-k[i] * t[i])
        else:                    # Ordre n
            base = 1.0 + (n - 1.0) * k[i] * t[i]
            if base < 0: base = 1e-10
            X_pred[i] = 1.0 - base**(1.0 / (1.0 - n))
            
    return X_pred

# ==========================================
# 2. CONFIGURATION DE L'ÉTUDE
# ==========================================
# C'est ici que tu choisis quel procédé tu étudies aujourd'hui !
PROCESSUS_CIBLE = "Neutral_Hydrolysis"  # Ex: "Acidic_Hydrolysis", "Alcoholysis"
CLE_JSON = "neutre"                     # Ex: "acide", "alcoolyse"

grille_temperatures_C = [220, 240, 260, 280]
grille_temps_h = [0.5, 1.0, 2.0, 4.0]
# On rajoute 0.0 dans les doses de catalyseur (très important pour l'hydrolyse neutre)
grille_catalyseur_wt = [0.0, 0.5, 2.0, 5.0] 

# ==========================================
# 3. ENTRAÎNEMENT DYNAMIQUE DE L'IA
# ==========================================
print(f"\n--- 1. Entraînement de l'IA pour {PROCESSUS_CIBLE} ---")
df, features, scaler, structure = charger_et_preparer_donnees()

if df is None:
    print("Erreur de chargement des données.")
    exit()

expert = entrainer_expert(PROCESSUS_CIBLE, df, features)
if expert is None:
    print(f"Échec de l'entraînement de l'IA pour {PROCESSUS_CIBLE}.")
    exit()

print("Modèle IA entraîné avec succès !")

# ==========================================
# 4. GÉNÉRATION DES DONNÉES PAR L'IA
# ==========================================
print("\n--- 2. Lancement du Laboratoire Virtuel (Génération des points) ---")

temps_list, temp_K_list, cat_list, rendements_predits = [], [], [], []

for T_C, t_h, cat_wt in itertools.product(grille_temperatures_C, grille_temps_h, grille_catalyseur_wt):
    
    # Création d'une ligne de données vide avec les bonnes colonnes pour tromper le Scaler
    ligne = {col: 0.0 for col in features}
    
    for col in features:
        col_low = col.lower()
        if 'temp' in col_low or 'ºc' in col_low: ligne[col] = T_C
        elif 'press' in col_low or 'kpa' in col_low: ligne[col] = 15000.0  # Pression arbitraire haute
        elif 'time' in col_low or 'h)' in col_low: ligne[col] = t_h
        elif 'liq' in col_low or 'ml/g' in col_low: ligne[col] = 5.0
        elif 'wt%' in col_low and 'add' not in col_low: ligne[col] = cat_wt
        
        # On active la colonne du procédé ciblé
        if PROCESSUS_CIBLE in col or PROCESSUS_CIBLE.replace(" ", "_") in col: 
            ligne[col] = 1.0

    # L'IA fait sa prédiction
    df_test = pd.DataFrame([ligne], columns=features)
    scaled = scaler.transform(df_test)
    rendement = expert.predict(pd.DataFrame(scaled, columns=features))[0]
    
    # Sécurité mathématique
    rendement = np.clip(rendement, 0.01, 0.999)
    
    temps_list.append(t_h)
    temp_K_list.append(T_C + 273.15)
    cat_list.append(cat_wt)
    rendements_predits.append(rendement)

print(f"{len(rendements_predits)} expériences virtuelles générées.")

# ==========================================
# 5. AJUSTEMENT CINÉTIQUE (SOLVEUR SCIPY)
# ==========================================
X_data = (np.array(temps_list), np.array(temp_K_list), np.array(cat_list))
Y_data = np.array(rendements_predits)

bornes = ([1e3, 40000, 0.1, 0.0], [1e12, 150000, 3.0, 2.0])
guess = [1e6, 80000, 1.0, 1.0]

print("\n--- 3. Optimisation Thermodynamique (Scipy) ---")
parametres_opti, _ = curve_fit(modele_cinetique_n_m, X_data, Y_data, p0=guess, bounds=bornes)

A_opt, Ea_opt, n_opt, m_opt = parametres_opti

print(f"Facteur pré-exponentiel (A) : {A_opt:.2e} h^-1")
print(f"Énergie d'activation (Ea)   : {Ea_opt/1000:.2f} kJ/mol")
print(f"Ordre de réaction (n)       : {n_opt:.2f}")
print(f"Ordre catalytique (m)       : {m_opt:.2f}")

# ==========================================
# 6. SAUVEGARDE INTELLIGENTE POUR BIOSTEAM
# ==========================================
fichier_json = os.path.join(DOSSIER_TRAVAIL, "parametres_cinetiques.json")

nouvelles_donnees = {
    "Ea_J_mol": float(Ea_opt),
    "A_h_1": float(A_opt),
    "n_ordre": float(n_opt),
    "m_catalyseur": float(m_opt)
}

# On lit l'existant s'il y en a un pour ne rien écraser
if os.path.exists(fichier_json):
    with open(fichier_json, "r") as f:
        base_de_donnees = json.load(f)
else:
    base_de_donnees = {}

# On met à jour uniquement le procédé ciblé
base_de_donnees[CLE_JSON] = nouvelles_donnees

# On sauvegarde tout
with open(fichier_json, "w") as f:
    json.dump(base_de_donnees, f, indent=4)

print(f"\nParamètres pour '{CLE_JSON}' fusionnés avec succès dans :")
print(f"-> {fichier_json}")