import pandas as pd
import os
import numpy as np
import json
import matplotlib.pyplot as plt
import math
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler
from pont_ia_biosteam import generer_flux_sortie_reacteur, orchestrateur_usine
import warnings
warnings.filterwarnings('ignore')

REACTEURS_BIBLIO = {
    "H3PO4": ["Closed loop reactor", "Steam-sparged stirred vessel"],
    "H3PW12O40 HPA": ["Batch"],
    "Fe- MOF-1": ["High-pressure reaction kettle"],
    "Fe- MOF-2": ["High-pressure reaction kettle"],
    "m-ZrO2": ["Autoclave batch"],
    "0": ["Batch", "Autoclave batch"], 
    "NaOH-KOH": ["Closed loop reactor"],
    "NaOH": ["Batch"],
    "HCl": ["Batch"],
    "H2SO4": ["Batch"],
    "MSA": ["Batch"],
    "p-TSA": ["Batch"],
    "BSA": ["Batch"]
}

# ---------------------------------------------------------------------
# SÉPARATION RIGOUREUSE DES SOLVANTS (THERMODYNAMIQUE)
# ---------------------------------------------------------------------
ALCOOLS_AUTOCLAVE = [
    "Methanol", "EtOH", "1-Propanol", "2-Propanol", 
    "1-Butanol", "2-Butanol", "2-Methyl-1-propanol", "2-Methyl-2-propanol"
]

ALCOOLS_BATCH = [
    "1-Pentanol", "2.2-Dimethyl-1-propanol", "1-Hexanol", "Cyclohexanol", 
    "1-Heptanol", "1-Octanol", "1-Decanol", "Dodecanol", "Cyclododecanol"
]
# ---------------------------------------------------------------------

DOSSIER_TRAVAIL = r"c:\Users\mathi\Desktop\Stage_Nylon\python\superstructura"

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
    
    # MODE ALCOOLYSE : Évaluation robuste OOB pour les petites bases de données
    if "alcoholysis" in proc_choisi.lower() or "glycolysis" in proc_choisi.lower():
        params = {'n_estimators': [100, 300], 'max_depth': [5, 7, 10], 'oob_score': [True]}
        rf = GridSearchCV(RandomForestRegressor(random_state=42), params, cv=3)
        rf.fit(X, y)
        print(f"Ready-to-use template : {rf.best_params_}")
        score = getattr(rf.best_estimator_, "oob_score_", 0.0)
        print(f" -> [Alcoholysis Mode] R² (Out-Of-Bag) : {score * 100:>6.2f} %")
        return rf.best_estimator_
    
    # MODE HYDROLYSE : Ancienne méthode sur 100% des données (Performances maximales)
    else:
        params = {'n_estimators': [100, 300], 'max_depth': [5, 7, 10], 'oob_score': [True]}
        rf = GridSearchCV(RandomForestRegressor(random_state=42), params, cv=3)
        rf.fit(X, y)
        print(f"Ready-to-use template : {rf.best_params_}")
        score_rf = getattr(rf.best_estimator_, "oob_score_", 0.0)
        print(f" -> [Hydrolysis Mode] R² (Out-Of-Bag) : {score_rf * 100:>6.2f} %")
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
        try:
            idx = int(input("Your choice (number) : "))
            proc_choisi = procs[idx]
        except (ValueError, IndexError):
            print(" Choix invalide. On prend le premier.")
            proc_choisi = procs[0]

        if proc_choisi not in experts:
            experts[proc_choisi] = entrainer_expert(proc_choisi, df, features)
            if experts[proc_choisi] is None: continue

        expert_actif = experts[proc_choisi]
        liste_cats = structure[proc_choisi]
        print(f"\n [BRANCHE] Step 2 : Select the catalyst :")
        for i, cat in enumerate(liste_cats): print(f" [{i}] {cat}")
        try:
            idx_cat = int(input("Your choice : "))
            cat_choisi = str(liste_cats[idx_cat])
        except:
            print("Invalid, defaulting to first.")
            cat_choisi = str(liste_cats[0])

        reacteurs_possibles = REACTEURS_BIBLIO.get(cat_choisi, ["Batch"])

        if len(reacteurs_possibles) == 1:
            reacteur_choisi = reacteurs_possibles[0]
            print(f" Reactor automatically assigned by the literature : {reacteur_choisi}")
        else:
            print(f"\n [ÉQUIPEMENT] literature offers several choice for '{cat_choisi}'.")
            for i, r in enumerate(reacteurs_possibles): 
                print(f"  [{i}] {r}")
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

        print("\n [CONDITIONS] Step 4: Enter the physical parameters:")

        if "H3PO4" in cat_choisi: t_min, t_max = 210.0, 370.0; p_min, p_max = 4000.0, 15000.0
        elif "ZrO2" in cat_choisi: t_min, t_max= 180.0, 230.0; p_min, p_max = 1000, 1000
        elif "HCl" in cat_choisi: t_min, t_max = 100.0, 105.0; p_min, p_max = 101.0, 101.1
        elif "Fe- MOF-1" in cat_choisi: t_min, t_max = 220,260; p_min, p_max = 2320, 4690
        elif "Fe- MOF-2" in cat_choisi: t_min, t_max = 220,260; p_min, p_max = 2320, 4690
        elif "H3PW12O40 HPA" in cat_choisi: t_min, t_max = 280,330; p_min, p_max = 9000, 9000
        elif "NaOH-KOH" in cat_choisi: t_min, t_max = 290,350; p_min, p_max = 101,101
        elif "NaOH" in cat_choisi: t_min, t_max = 225,270; p_min, p_max = 0.4,5000
        else: t_min, t_max = 180, 620; p_min, p_max = 103, 64200

        while True:
            temp = float(input(f" -> Temperature (°C) [Strict limits: {t_min:.1f} à {t_max:.1f}] : "))
            if t_min <= temp <= t_max: break
            print("  Error: Temperature outside the range specified in the literature.")

        p_sat_bar = math.exp(11.68 - 3816.44 / (temp + 227.02))
        p_sat_kpa = p_sat_bar * 100
        
        if "Acidic" in str(proc_choisi) or "H3PO4" in str(cat_choisi):
            print(f" -> [INFO] Vapor Stripping detected. Water boils at {p_sat_kpa:.1f} KPa at {temp}°C.")
            while True:
                pression = float(input(f" -> Pressure (KPa) [Must be < {p_sat_kpa:.1f} KPa to allow boiling] : "))
                if pression < p_sat_kpa: break
                print(f"  Error: Pressure too high ! Water will not boil.")
        else:
            pression = p_sat_kpa + 2000.0
            print(f" -> [AUTO] Pressure optimized to : {pression:.1f} KPa (to maintain stable liquid phase)")

        time = float(input(" -> reaction time (h) [Min: 0.1 | Max: 72.0] : "))
        
        if cat_choisi == "0":
            cat_val = 0.0
            print(" -> Amount of catalyst (wt%) : 0 (No catalyst selected)")
        else:
            cat_val = float(input(f" -> Amount of catalyst {cat_choisi} (wt%) [Min: 0 | Max: 100] : "))
            
        liq_ratio = float(input(" -> liquid ratio (mL/g) [Min: 0.1 | Max: 200] : "))
        
        try:
            purete_pa6 = float(input(" -> plastic waste purity (%) [Min: 10 | Max: 100] : "))
        except ValueError:
            purete_pa6 = 60.0
        
        # ---------------------------------------------------------------------
        # LE MENU S'ADAPTE AU RÉACTEUR
        # ---------------------------------------------------------------------
        if "alcoholysis" in str(proc_choisi).lower() or "glycolysis" in str(proc_choisi).lower():
            print(f"\n [SOLVENT] Vous utilisez un {reacteur_choisi}. Seuls les solvants compatibles sont affichés.")
            
            liste_a_afficher = ALCOOLS_AUTOCLAVE if "Autoclave" in reacteur_choisi else ALCOOLS_BATCH
            
            for i, alc in enumerate(liste_a_afficher): 
                print(f"  [{i}] {alc}")
                
            while True:
                try:
                    idx_alc = int(input(" -> Choisissez votre alcool (numéro) : "))
                    if 0 <= idx_alc < len(liste_a_afficher):
                        nom_solvant = liste_a_afficher[idx_alc]
                        print(f" -> Solvant sélectionné : {nom_solvant}")
                        break
                    else: print("  Numéro invalide.")
                except ValueError: print("  Erreur : Veuillez entrer un numéro.")
        else:
            nom_solvant = 'Water'

        nom_additif = input(" -> Additive name (Type '0' if no additive, e.g., 'CL') : ")
        if nom_additif != '0':
            additif_val = float(input(f" -> Amount of additive {nom_additif} (wt%) : "))
        else:
            additif_val = 0.0

        # Remplissage sécurisé de la ligne IA avec str() partout
        prototype_ligne = {col: 0.0 for col in features}
        
        for col in features:
            col_low = col.lower()
            if 'reactor' not in col_low and 'process' not in col_low and 'type' not in col_low:
                if 'temp' in col_low or 'ºc' in col_low: prototype_ligne[col] = temp
                elif 'press' in col_low or 'kpa' in col_low: prototype_ligne[col] = pression
                elif 'time' in col_low or 'h)' in col_low: prototype_ligne[col] = time
                elif 'liq' in col_low or 'ml/g' in col_low: prototype_ligne[col] = liq_ratio
                elif 'wt%' in col_low and 'add' not in col_low: prototype_ligne[col] = cat_val
                elif 'wt%' in col_low and 'add' in col_low: prototype_ligne[col] = additif_val

            if str(proc_choisi) in col or str(proc_choisi).replace(" ", "_") in col: prototype_ligne[col] = 1.0
            if str(cat_choisi) in col and 'type' in col_low: prototype_ligne[col] = 1.0
            if str(reacteur_choisi) in col: prototype_ligne[col] = 1.0
            if nom_additif != '0' and col.endswith(str(nom_additif)): prototype_ligne[col] = 1.0
            if nom_solvant != 'Water' and str(nom_solvant).lower() in col_low: prototype_ligne[col] = 1.0

            if 'cl' in col_low and not 'hcl' in col_low and not 'auto' in col_low and not 'closed' in col_low:
                if 'target' in col_low or 'product' in col_low: prototype_ligne[col] = 1.0
            if 'open' in col_low or 'chen' in col_low or 'src_chen' in col_low:
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

        if input("\n Run the BioSTEAM plant simulation for this prediction ? (o/n) : ").lower() == 'o':
            print(" LAUNCH OF THE INDUSTRIAL DIGITAL TWIN")
            masse_dechets_entrants = 1000.0 
            masse_nylon_pur = masse_dechets_entrants * (purete_pa6 / 100.0)
            masse_fibres = masse_dechets_entrants - masse_nylon_pur
            masse_eau_usine = masse_nylon_pur * liq_ratio
            masse_cat_usine = masse_nylon_pur * (cat_val / 100.0)

            rendement_fraction = pred_rf / 100.0
            conversion_estimee = min(1.0, rendement_fraction + 0.05) 
            predictions_ia = {'Conversion': conversion_estimee, 'Rendement_CPL': rendement_fraction}

            flux_usine = generer_flux_sortie_reacteur(
                masse_pa6_in=masse_nylon_pur, masse_eau_in=masse_eau_usine,
                type_catalyseur=cat_choisi, masse_catalyseur=masse_cat_usine,
                predictions_IA=predictions_ia, masse_fibres=masse_fibres, nom_solvant=nom_solvant
            )

            try:
                kpi_usine = orchestrateur_usine(
                    flux_biosteam=flux_usine, type_catalyseur=cat_choisi, type_reacteur=reacteur_choisi,
                    temp_user=temp, pression_user=pression, temps_sejour=time, nom_solvant=nom_solvant
                )
                
                revenu_cpl = flux_usine.get('Caprolactam', 0) * 2.0 
                cout_vapeur = kpi_usine.get('Cout_Vapeur_hr', 0.0)
                cout_eau = kpi_usine.get('Cout_Eau_Froide_hr', 0.0)
                cout_solvant = kpi_usine.get('Cout_Solvant_hr', 0.0)
                cout_equipement = kpi_usine.get('Cout_Equipement_hr', 0.0)
                cout_reactifs = kpi_usine.get('Cout_Reactifs_hr', 0.0) 
                cout_catalyseur = kpi_usine.get('Cout_Catalyseur_hr', 0.0) 
                
                if kpi_usine.get('Type_Usine', "") == "Neutre":
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

                couts_totaux = cout_vapeur + cout_eau + cout_solvant + cout_equipement + cout_reactifs + cout_catalyseur + cout_dechets
                profit = revenu_cpl - couts_totaux
                
                if flux_usine.get('Caprolactam', 0) > 0:
                    msp = couts_totaux / flux_usine.get('Caprolactam', 0)
                else:
                    msp = float('inf')

                print(f" -> Active Architecture       : Usine {kpi_usine.get('Type_Usine', 'Inconnu')}")
                print(f" -> CPL produced              : {flux_usine.get('Caprolactam', 0):.1f} kg/h")
                print(f" -> 💰 MSP (Min Selling Price) : {msp:.2f} €/kg")
                print(f" -> Purity of Caprolactam     : {kpi_usine.get('Purete_CPL', 0):.2f} %")
                print(f" -> {nom_dechet_1}{nom_dechet_2}")
                print(f" -> CAPEX                     : {cout_equipement:.2f} € / heure")
                print(f" -> OPEX Hot Steam            : {cout_vapeur:.2f} € / heure")
                print(f" -> OPEX Cooling Water        : {cout_eau:.2f} € / heure")
                print(f" -> OPEX Solvent              : {cout_solvant:.2f} € / heure")
                if cout_reactifs > 0: print(f" -> OPEX Neutralizing Reagents: {cout_reactifs:.2f} € / heure")
                if cout_catalyseur > 0: print(f" -> OPEX Purchase catalyst      : {cout_catalyseur:.2f} € / heure")
                print(f" -> ESTIMATED HOURLY PROFIT      : {profit:.2f} € / heure")
                
            except Exception as e:
                print(f"Error while launching the factory : {e}")

        if input("\n Plot the trend over time? (o/n) : ").lower() == 'o':
            print(" -> Generation of virtual predictions (0 à 72h)...")
            temps_virtuels = np.linspace(0.0, 72.0, 100)
            preds_rf_plot = []
            for t in temps_virtuels:
                if t == 0.0: preds_rf_plot.append(0.0) 
                else:
                    ligne_virtuelle = prototype_ligne.copy()
                    for col in features:
                        if 'time' in col.lower() or 'h)' in col.lower():
                            if 'reactor' not in col.lower() and 'process' not in col.lower():
                                ligne_virtuelle[col] = t
                    df_virtuel = pd.DataFrame([ligne_virtuelle], columns=features)
                    scaled_virtuel = scaler.transform(df_virtuel)
                    df_virtuel_scaled = pd.DataFrame(scaled_virtuel, columns=features)
                    prediction = expert_actif.predict(df_virtuel_scaled)[0]
                    preds_rf_plot.append(max(0.0, prediction))

            plt.figure(figsize=(10, 6))
            plt.plot(temps_virtuels, preds_rf_plot, label='Random Forest', color='orange', linestyle='--', linewidth=2)
            plt.scatter([time], [pred_rf], color='orange', zorder=5, s=100, label='Point actuel')
            titre_alcool = f" | Solvent: {nom_solvant}" if nom_solvant != 'Water' else ""
            plt.title(f"Simulateur de Cinétique : {str(proc_choisi).replace('_', ' ')}\n(Temp: {temp}°C | Pression: {pression} KPa | Cat: {cat_val}wt%{titre_alcool})", fontsize=12)
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