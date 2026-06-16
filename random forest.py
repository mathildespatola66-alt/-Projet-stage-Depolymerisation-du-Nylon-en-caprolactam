import pandas as pd
import os
import numpy as np
import json
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from pont_ia_biosteam import generer_flux_sortie_reacteur, orchestrateur_usine

# =========================================================================
# DICTIONNAIRE BIBLIOGRAPHIQUE RÉEL
# =========================================================================
REACTEURS_BIBLIO = {
    "H3PO4": ["Closed loop reactor", "Steam-sparged stirred vessel"],
    "H3PW12O40 HPA": ["Batch"],
    
    "Fe- MOF-1": ["High-pressure reaction kettle"],
    "Fe- MOF-2": ["High-pressure reaction kettle"],
    "m-ZrO2": ["Autoclave batch"],
    "0": ["Bench-top stirred reactor", "Batch", "Autoclave batch"],
    
    "NaOH-KOH": ["Closed loop reactor"],
    "NaOH": ["Batch"]
}

DOSSIER_TRAVAIL = r"c:\Users\mathi\Desktop\Stage_Nylon\python\process"

def charger_et_preparer_donnees():
    print("\n--- READING AND SORTING DATA ---")
    chemin_data = os.path.join(DOSSIER_TRAVAIL, 'donnees_preparees.csv')
    chemin_json = os.path.join(DOSSIER_TRAVAIL, "structure_arbre_chimique.json")
    
    if not os.path.exists(chemin_data) or not os.path.exists(chemin_json):
        print("Files not found. Check the paths.")
        return None, None, None, None

    df = pd.read_csv(chemin_data, sep=';')
    with open(chemin_json, "r", encoding="utf-8") as f:
        structure = json.load(f)
    
    features = [c for c in df.columns if c != 'Yield (%)']
    scaler = MinMaxScaler()
    df[features] = scaler.fit_transform(df[features])
    
    return df, features, scaler, structure

def entrainer_expert(proc_choisi, df, features):
    print(f"\nTraining a Random Forest for : {proc_choisi.replace('_', ' ')}")
    col_proc = f"Process_Type_{proc_choisi}"
    df_expert = df[df[col_proc] == 1]
    
    if len(df_expert) < 5:
        print("Not enough data for this process.")
        return None

    X = df_expert[features]
    y = df_expert['Yield (%)']
    
    params = {'n_estimators': [100, 300], 'max_depth': [5, 7, 10]}
    rf = GridSearchCV(RandomForestRegressor(random_state=42), params, cv=3)
    rf.fit(X, y)
    print(f"Ready-to-use template : {rf.best_params_}")
    return rf.best_estimator_

def lancer_assistant():
    df, features, scaler, structure = charger_et_preparer_donnees()
    if df is None: return
    
    experts = {}

    while True:
        print("\n" + "="*60)
        print(" CHEMICAL ASSISTANT - RANDOM FOREST ")
        print("="*60)

        procs = list(structure.keys())
        print("\n Step 1: Which process do you want to study?")
        for i, p in enumerate(procs): print(f" [{i}] {p.replace('_', ' ')}")
        idx = int(input("Your choice (number) : "))
        proc_choisi = procs[idx]

        if proc_choisi not in experts:
            experts[proc_choisi] = entrainer_expert(proc_choisi, df, features)
            if experts[proc_choisi] is None: continue

        expert_actif = experts[proc_choisi]

        liste_cats = structure[proc_choisi]
        print(f"\n [BRANCHE] Step 2 : Select the catalyst :")
        for i, cat in enumerate(liste_cats): print(f" [{i}] {cat}")
        idx_cat = int(input("Your choice : "))
        cat_choisi = liste_cats[idx_cat]

        reacteurs_possibles = REACTEURS_BIBLIO.get(cat_choisi, ["Batch"])

        if len(reacteurs_possibles) == 1:
            reacteur_choisi = reacteurs_possibles[0]
            print(f" Reactor automatically assigned by the literature : {reacteur_choisi}")
        else:
            print(f"\n [ÉQUIPEMENT] literature offers several choice for '{cat_choisi}'.")
            for i, r in enumerate(reacteurs_possibles): 
                print(f"  [{i}] {r}")
            
            # Boucle de sécurité pour éviter les fautes de frappe
            while True:
                try:
                    idx_r = int(input("  Which one would you like to simulate? (number): "))
                    if 0 <= idx_r < len(reacteurs_possibles):
                        reacteur_choisi = reacteurs_possibles[idx_r]
                        print(f" Selected reactor: {reacteur_choisi}")
                        break
                    else:
                        print("   Invalid number. Please select a number from the list.")
                except ValueError:
                    print("   Error: Please enter a NUMBER (such as 0 or 1), not a letter!")


        # =====================================================================
        # [CONDITIONS] Étape 4 : Entrées physiques et Limites Littérature
        # =====================================================================
        print("\n [CONDITIONS] Step 4: Enter the physical parameters:")


        if "H3PO4" in str(cat_choisi):
            t_min, t_max = 210.0, 370.0
            # CORRECTION : On passe les limites en KPa pour empêcher l'eau de bouillir !
            # (1000 PSI = ~6900 KPa, on met une limite haute de sécurité à 15000 KPa)
            p_min, p_max = 4000.0, 15000.0
        elif "ZrO2" in str(cat_choisi):
            t_min, t_max= 180.0, 230.0
            p_min, p_max = 1000, 1000
        elif "HCl" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 100.0, 105.0
            p_min, p_max = 101.0, 101.1
        elif "Fe- MOF-1" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 220,260
            p_min, p_max = 2320, 4690
        elif "Fe- MOF-2" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 220,260
            p_min, p_max = 2320, 4690
        elif "H3PW12O40 HPA" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 280,330
            p_min, p_max = 9000, 9000
        elif "NaOH-KOH" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 290,350
            p_min, p_max = 101,101
        elif "NaOH" in str(cat_choisi) in str(cat_choisi):
            t_min, t_max = 225,270
            p_min, p_max = 0.4,5000
        else:
            t_min, t_max = 180, 620
            p_min, p_max = 103, 64200

     

        # 1. Demande de la Température
        while True:
            temp = float(input(f" -> Temperature (°C) [Strict limits: {t_min:.1f} à {t_max:.1f}] : "))
            if t_min <= temp <= t_max: break
            print("  Error: Temperature outside the range specified in the literature.")

        # =====================================================================
        # 2. CALCUL AUTOMATIQUE DE LA PRESSION (Loi d'Antoine Thermodynamique)
        # =====================================================================
        import math
        
        # Calcul de la Pression de Saturation de l'eau (en bar)
        p_sat_bar = math.exp(11.68 - 3816.44 / (temp + 227.02))
        
        if proc_choisi == "Acidic Hydrolysis":
            # Le procédé acide a BESOIN de bouillir pour extraire le monomère (Stripping)
            # On force une basse pression (ex: 400 KPa = 4 bar)
            pression = 400.0
            print(f" -> [AUTO] Pressure optimized to : {pression:.1f} KPa (Vapor Stripping enabled)")
        else:
            # Procédés Neutre et Basique : On force l'état liquide avec une marge de 20 bar
            pression = (p_sat_bar * 100) + 2000.0
            print(f" -> [AUTO] Pressure optimized to : {pression:.1f} KPa (to maintain stable liquid phase)")
        # =====================================================================
        # =====================================================================

        
        time = float(input(" -> reaction time (h) [Min: 0.1 | Max: 72.0] : "))
        cat_val = float(input(" -> amount of catalyst (wt%) [Min: 0 | Max: 100] : "))
        liq_ratio = float(input(" -> liquid ratio (mL/g) [Min: 0.1 | Max: 200] : "))
        purete_pa6 = float(input(" -> plastic waste purity (%) [Min: 10 | Max: 100] : "))

        
        prototype_ligne = {col: 0.0 for col in features}
        
        for col in features:
            if 'Temperature' in col or 'ºC' in col: prototype_ligne[col] = temp
            elif 'Pressure' in col or 'KPa' in col: prototype_ligne[col] = pression
            elif 'Time' in col or '(h)' in col: prototype_ligne[col] = time
            elif 'wt%' in col: prototype_ligne[col] = cat_val
            elif 'liq' in col.lower() or 'ml/g' in col.lower(): prototype_ligne[col] = liq_ratio

            if str(proc_choisi).lower() in col.lower() or str(cat_choisi).lower() in col.lower():
                prototype_ligne[col] = 1.0
            if str(reacteur_choisi).lower() in col.lower():
                prototype_ligne[col] = 1.0
            

        if time == 0.0:
            pred_rf = 0.0
        else:
            df_custom = pd.DataFrame([prototype_ligne], columns=features)
            custom_scaled = scaler.transform(df_custom)
            df_custom_scaled = pd.DataFrame(custom_scaled, columns=features)
            pred_rf = expert_actif.predict(df_custom_scaled)[0]
            pred_rf = max(0.0, pred_rf)

        print("\n" + "="*60)
        print(f"OPINION OF THE AI EXPERT IN : {str(proc_choisi).replace('_', ' ').upper()}")
        print(f"  [Rendement Prédit par Random Forest] : {pred_rf:.2f} %")
        print("="*60)

        # =====================================================================
        # COUPLAGE DIRECT : IA -> BIOSTEAM
        # =====================================================================
        if input("\n Run the BioSTEAM plant simulation for this prediction ? (o/n) : ").lower() == 'o':
            print(" LAUNCH OF THE INDUSTRIAL DIGITAL TWIN")
            
            masse_dechets_entrants = 1000.0 
            masse_nylon_pur = masse_dechets_entrants * (purete_pa6 / 100.0)
            masse_fibres = masse_dechets_entrants - masse_nylon_pur
            
            masse_eau_usine = masse_nylon_pur * liq_ratio
            masse_cat_usine = masse_nylon_pur * (cat_val / 100.0)

            rendement_fraction = pred_rf / 100.0
            conversion_estimee = min(1.0, rendement_fraction + 0.05) 

            predictions_ia = {
                'Conversion': conversion_estimee,
                'Rendement_CPL': rendement_fraction
            }

            flux_usine = generer_flux_sortie_reacteur(
                masse_pa6_in=masse_nylon_pur,
                masse_eau_in=masse_eau_usine,
                type_catalyseur=cat_choisi,
                masse_catalyseur=masse_cat_usine,
                predictions_IA=predictions_ia,
                masse_fibres=masse_fibres 
            )

            try:
                # N'oublie pas que orchestrateur_usine doit aussi accepter et relayer 'temps_sejour' à simuler_usine_neutre
                kpi_usine = orchestrateur_usine(
                    flux_biosteam=flux_usine, 
                    type_catalyseur=cat_choisi, 
                    type_reacteur=reacteur_choisi,
                    temp_user=temp,          
                    pression_user=pression,
                    temps_sejour=time  
                )
                
                
                revenu_cpl = flux_usine['Caprolactam'] * 2.0 
                
                cout_vapeur = kpi_usine.get('Cout_Vapeur_hr', 0.0)
                cout_eau = kpi_usine.get('Cout_Eau_Froide_hr', 0.0)
                cout_solvant = kpi_usine.get('Cout_Solvant_hr', 0.0)
                cout_equipement = kpi_usine.get('Cout_Equipement_hr', 0.0)
                cout_reactifs = kpi_usine.get('Cout_Reactifs_hr', 0.0) 
                cout_catalyseur = kpi_usine.get('Cout_Catalyseur_hr', 0.0) 
                
                if kpi_usine['Type_Usine'] == "Neutre":
                    eau_usee = kpi_usine.get('Dechets_Totaux_kg', 0)
                    cout_dechets = eau_usee * 0.002 
                    nom_dechet_1 = f"Eaux usées + Fibres        : {eau_usee:.1f} kg/h"
                    nom_dechet_2 = ""
                else:
                    eau_usee = kpi_usine.get('Eaux_Usees_kg', 0)
                    sels = kpi_usine.get('Dechets_Toxiques_kg', 0)
                    cout_dechets = (eau_usee * 0.002) + (sels * 0.50) 
                    nom_dechet_1 = f"Eaux usées classiques      : {eau_usee:.1f} kg/h"
                    nom_dechet_2 = f"\n -> Sels + Déchets toxiques    : {sels:.1f} kg/h"

                
                profit = revenu_cpl - cout_dechets - cout_vapeur - cout_eau - cout_equipement - cout_solvant - cout_reactifs - cout_catalyseur

                
                print(f" -> Active Architecture       : Usine {kpi_usine['Type_Usine']}")
                print(f" -> CPL produced                : {flux_usine['Caprolactam']:.1f} kg/h (depuis {masse_nylon_pur:.0f} kg de PA6)")
                print(f" -> Purity of Caprolactam     : {kpi_usine['Purete_CPL']:.2f} %")
                print(f" -> {nom_dechet_1}{nom_dechet_2}")
                print(f" -> CAPEX       : {cout_equipement:.2f} € / heure")
                print(f" -> OPEX Hot Steam         : {cout_vapeur:.2f} € / heure")
                print(f" -> OPEX Cooling Water: {cout_eau:.2f} € / heure")
                print(f" -> OPEX Solvent    : {cout_solvant:.2f} € / heure")
                
                if cout_reactifs > 0:
                    print(f" -> OPEX Neutralizing Reagents: {cout_reactifs:.2f} € / heure")
                
                if cout_catalyseur > 0:
                    print(f" -> OPEX Purchase catalyst      : {cout_catalyseur:.2f} € / heure")
                    
                print(f" -> ESTIMATED HOURLY PROFIT      : {profit:.2f} € / heure")
                
            except Exception as e:
                print(f"Error while launching the factory : {e}")

        # =====================================================================
        # 7. TRACÉ DU GRAPHIQUE
        # =====================================================================
        if input("\n Plot the trend over time? (o/n) : ").lower() == 'o':
            print(" -> Generation of virtual predictions (0 à 72h)...")
            
            temps_virtuels = np.linspace(0.0, 72.0, 100)
            preds_rf_plot = []
            
            for t in temps_virtuels:
                if t == 0.0:
                    preds_rf_plot.append(0.0) 
                else:
                    ligne_virtuelle = prototype_ligne.copy()
                    for col in features:
                        if 'time' in col.lower() or 'h)' in col.lower():
                            ligne_virtuelle[col] = t
                    
                    df_virtuel = pd.DataFrame([ligne_virtuelle], columns=features)
                    scaled_virtuel = scaler.transform(df_virtuel)
                    df_virtuel_scaled = pd.DataFrame(scaled_virtuel, columns=features)
                    
                    prediction = expert_actif.predict(df_virtuel_scaled)[0]
                    preds_rf_plot.append(max(0.0, prediction))

            plt.figure(figsize=(10, 6))
            plt.plot(temps_virtuels, preds_rf_plot, label='Random Forest', color='orange', linestyle='--', linewidth=2)
            plt.scatter([time], [pred_rf], color='orange', zorder=5, s=100, label='Point actuel')
            
            plt.title(f"Simulateur de Cinétique : {str(proc_choisi).replace('_', ' ')}\n(Temp: {temp}°C | Pression: {pression} KPa | Cat: {cat_val}wt%)", fontsize=12)
            plt.xlabel("Temps de réaction (Heures)", fontsize=11)
            plt.ylabel("Rendement Prédit (%)", fontsize=11)
            plt.ylim(-5, 105)
            plt.grid(True, linestyle=':', alpha=0.7)
            plt.legend(fontsize=10)
            plt.show()

        if input("\n Try a different configuration? (o/n) : ").lower() != 'o':
            break

if __name__ == "__main__":
    lancer_assistant()