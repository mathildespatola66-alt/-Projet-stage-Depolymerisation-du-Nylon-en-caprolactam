import numpy as np
import pandas as pd
import random
import json
import os
from scipy.optimize import least_squares

# Importation depuis tes propres scripts
from random_forest import charger_et_preparer_donnees, entrainer_expert, DOSSIER_TRAVAIL
from superstructure_globale import creer_ligne_ia, obtenir_limites_physiques, REACTEURS_BIBLIO

# Liste des catalyseurs qui n'ont pas de cinétique dans ton JSON
CATALYSEURS_MANQUANTS = {
    "Acidic_Hydrolysis": ["HCl", "H3PO4"],
    "Alkaline_Hydrolysis": ["NaOH-KOH"],
    "Neutral_Hydrolysis":["m-ZrO2"]
}

POINTS_VIRTUELS_PAR_CAT = 500

def ajuster_cinetique_sur_ia(T_c, t_h, cat_wt, X_rf):
    """Trouve les constantes (A, Ea, n, m) qui collent le mieux aux prédictions de l'IA."""
    T_k = T_c + 273.15
    
    def residus(params):
        lnA, Ea, n, m = params
        A = np.exp(lnA)
        
        # Loi d'Arrhenius
        k = A * np.exp(-Ea / (8.314 * T_k)) * (np.clip(cat_wt, 1e-5, None)**m)
        
        # Intégration de l'ordre n (identique au script du professeur)
        eps = 1.0 - n
        eps_sur = np.where(np.abs(eps) < 1e-10, 1e-10, eps)
        arg = np.clip(-eps_sur * k * t_h, -1.0 + 1e-10, None)
        C = np.exp(np.log1p(arg) / eps_sur)
        X_pred = 1.0 - C
        
        return X_pred - X_rf

    # [ln(A), Ea (J/mol), n, m]
    x0 = [np.log(1e9), 80000.0, 1.5, 0.5]
    # Limites pour forcer le solveur à trouver des valeurs physiquement possibles
    limites = ([0.0, 10000.0, 0.1, -1.0], [50.0, 300000.0, 4.0, 3.0])
    
    res = least_squares(residus, x0, bounds=limites, method='trf')
    return res.x

def generer_et_ajuster():
    print("--- DÉMARRAGE DU GÉNÉRATEUR DE CINÉTIQUE PAR IA ---")
    df, features, scaler, structure = charger_et_preparer_donnees()
    
    resultats_json = {}
    
    for proc, catalyseurs in CATALYSEURS_MANQUANTS.items():
        expert = entrainer_expert(proc, df, features)
        if expert is None:
            continue
            
        for cat in catalyseurs:
            print(f"\n[Génération] Création de {POINTS_VIRTUELS_PAR_CAT} points virtuels pour {cat} ({proc})...")
            
            T_list, t_list, cat_wt_list, X_rf_list = [], [], [], []
            reacteur = REACTEURS_BIBLIO.get(cat, ["Batch"])[0]
            t_min, t_max, p_min, p_max, _ = obtenir_limites_physiques(proc, cat, reacteur)
            
            # 1. Génération des données avec la Forêt Aléatoire
            lignes_df = []
            for _ in range(POINTS_VIRTUELS_PAR_CAT):
                temp = random.uniform(t_min, t_max)
                press = random.uniform(p_min, p_max)
                time = random.uniform(0.5, 6.0)
                cat_val = random.uniform(0.5, 8.0)
                liq_ratio = random.uniform(1.0, 5.0)
                
                T_list.append(temp)
                t_list.append(time)
                cat_wt_list.append(cat_val)
                
                ligne = creer_ligne_ia(features, proc, cat, reacteur, temp, press, time, cat_val, "Water", "0", 0.0, liq_ratio)
                lignes_df.append(ligne)
                
            df_test = pd.DataFrame(lignes_df, columns=features)
            scaled = scaler.transform(df_test)
            rendements_predits = expert.predict(pd.DataFrame(scaled, columns=features))
            
            # Conversion en fraction (0.0 à 1.0) et filtrage des valeurs extrêmes pour le solveur
            X_rf_list = np.clip(rendements_predits / 100.0, 0.01, 0.99)
            
            # 2. Ajustement mathématique (Solveur)
            print(f"[{cat}] Ajustement de la loi de vitesse (Arrhenius) en cours...")
            parametres_optimises = ajuster_cinetique_sur_ia(
                np.array(T_list), np.array(t_list), np.array(cat_wt_list), X_rf_list
            )
            
            lnA, Ea, n, m = parametres_optimises
            A = float(np.exp(lnA))
            
            if "Acidic" in proc:
                cle_json = f"acide_{cat.lower()}"
            elif "Neutral" in proc:
                cle_json = f"neutre_{cat.lower()}"
            else:
                cle_json = f"alcaline_{cat.lower()}"
            resultats_json[cle_json] = {
                "Ea_J_mol": round(float(Ea), 1),
                "A_h_1": float(f"{A:.3e}"),
                "n_ordre": round(float(n), 3),
                "m_catalyseur": round(float(m), 3)
            }
            
            print(f"  -> A = {A:.2e} h-1 | Ea = {Ea/1000:.1f} kJ/mol | n = {n:.2f} | m = {m:.2f}")

    # 3. Sauvegarde
    chemin_sortie = os.path.join(DOSSIER_TRAVAIL, 'constantes_ia_generees.json')
    with open(chemin_sortie, 'w', encoding='utf-8') as f:
        json.dump(resultats_json, f, indent=4)
        
    print(f"\n Terminé ! Les constantes ont été sauvegardées dans : {chemin_sortie}")
    print("Tu peux maintenant ouvrir ce fichier et copier les blocs dans ton 'parametres_cinetiques.json' principal.")

if __name__ == "__main__":
    generer_et_ajuster()