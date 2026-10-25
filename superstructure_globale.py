import pandas as pd
import numpy as np
import random
import math
import sys
import os
import json
import warnings
import matplotlib.pyplot as plt

print("--- Démarrage du script Superstructure ---")

from pont_ia_biosteam import obtenir_params_pour_reacteur, orchestrateur_usine
from random_forest import charger_et_preparer_donnees, entrainer_expert, DOSSIER_TRAVAIL
print("--- Imports réussis ---")

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
    "MSA": ["Batch"],
    "p-TSA": ["Batch"],
    "BSA": ["Batch"]
}

ALCOOLS_AUTOCLAVE = [
    "Methanol", "EtOH", "1-Propanol", "2-Propanol", 
    "1-Butanol", "2-Butanol", "2-Methyl-1-propanol", "2-Methyl-2-propanol"
]

ALCOOLS_BATCH = [
    "1-Pentanol", "2.2-Dimethyl-1-propanol", "1-Hexanol", "Cyclohexanol", 
    "1-Heptanol", "1-Octanol", "1-Decanol", "Dodecanol", "Cyclododecanol"
]

# =========================================================================
# LES LOIS DE LA PHYSIQUE ET DE LA MÉCANIQUE
# =========================================================================
def obtenir_limites_physiques(proc, cat, reacteur):
    is_liquid_water = False
    p_min, p_max = 101.3, 101.3
    t_min, t_max = 200.0, 300.0

    if "alcoholysis" in proc.lower() or "glycolysis" in proc.lower():
        if "Autoclave" in reacteur:
            t_min, t_max = 250.0, 400.0
            p_min, p_max = 5000.0, 30000.0 
        else:
            t_min, t_max = 150.0, 260.0 
            p_min, p_max = 101.3, 500.0 
            
    elif "Acidic" in proc:
        if "H3PO4" in cat:
            t_min, t_max = 210.0, 370.0
            if "Steam-sparged" in reacteur:
                p_min, p_max = 500.0, 3500.0 
            else:
                is_liquid_water = True 
        elif cat in ["HCl", "H2SO4", "MSA", "p-TSA", "BSA"]:
            t_min, t_max = 100.0, 150.0
            p_min, p_max = 101.3, 500.0
        elif "HPA" in cat:
            t_min, t_max = 280.0, 330.0
            p_min, p_max = 8000.0, 10000.0
            
    elif "Neutral" in proc:
        if cat == "0":
            t_min, t_max = 250.0, 400.0
            is_liquid_water = True
        elif "ZrO2" in cat:
            t_min, t_max = 180.0, 230.0
            is_liquid_water = True
        elif "MOF" in cat:
            t_min, t_max = 220.0, 260.0
            p_min, p_max = 2000.0, 5000.0
            
    elif "Alkaline" in proc:
        is_liquid_water = True
        if "NaOH-KOH" in cat:
            t_min, t_max = 290.0, 350.0
        else:
            t_min, t_max = 225.0, 270.0
            
    return t_min, t_max, p_min, p_max, is_liquid_water

def creer_ligne_ia(features, proc, cat, reacteur, temp, press, time, cat_val, solvant, additif_nom, additif_val, liq_ratio):
    ligne = {col: 0.0 for col in features}
    for col in features:
        col_low = col.lower()
        if 'reactor' not in col_low and 'process' not in col_low and 'type' not in col_low:
            if 'temp' in col_low or 'ºc' in col_low: ligne[col] = temp
            elif 'press' in col_low or 'kpa' in col_low: ligne[col] = press
            elif 'time' in col_low or 'h)' in col_low: ligne[col] = time
            elif 'liq' in col_low or 'ml/g' in col_low: ligne[col] = liq_ratio
            elif 'wt%' in col_low and 'add' not in col_low: ligne[col] = cat_val
            elif 'wt%' in col_low and 'add' in col_low: ligne[col] = additif_val

        if str(proc) in col or str(proc).replace(" ", "_") in col: ligne[col] = 1.0
        if str(cat) in col and 'type' in col_low: ligne[col] = 1.0
        if str(reacteur) in col: ligne[col] = 1.0
        if str(solvant) != 'Water' and str(solvant).lower() in col_low: ligne[col] = 1.0
        if str(additif_nom) != '0' and col.endswith(str(additif_nom)): ligne[col] = 1.0
        
        if 'cl' in col_low and not 'hcl' in col_low and not 'auto' in col_low and not 'closed' in col_low:
            if 'target' in col_low or 'product' in col_low: ligne[col] = 1.0
    return ligne

def tracer_front_pareto(resultats):
    print("\n 📊 Génération du diagramme de Pareto (MSP vs GWP)...")
    if not resultats:
        print("❌ Aucun résultat à afficher.")
        return

    df = pd.DataFrame(resultats)
    
    if "MSP (€/kg)" in df.columns and "GWP_kgCO2_per_kgCPL" in df.columns:
        df = df[df["MSP (€/kg)"] != float('inf')]
        df = df[df["GWP_kgCO2_per_kgCPL"] != float('inf')]
    else:
        print("❌ Données incomplètes pour tracer le Pareto.")
        return

    if df.empty:
        print("❌ Aucune donnée valide pour le Pareto.")
        return

    df_sorted = df.sort_values(by="MSP (€/kg)")
    pareto_front = []
    min_gwp = float('inf')

    for index, row in df_sorted.iterrows():
        if row["GWP_kgCO2_per_kgCPL"] < min_gwp:
            pareto_front.append(row)
            min_gwp = row["GWP_kgCO2_per_kgCPL"]

    pareto_df = pd.DataFrame(pareto_front)

    plt.figure(figsize=(10, 6))
    plt.scatter(df["MSP (€/kg)"], df["GWP_kgCO2_per_kgCPL"], color='#94a3b8', alpha=0.6, label='Scénarios explorés', s=50)
    plt.plot(pareto_df["MSP (€/kg)"], pareto_df["GWP_kgCO2_per_kgCPL"], color='#ef4444', marker='o', linestyle='-', linewidth=2, label='Front de Pareto (Optimums)', markersize=8)

    for i, row in pareto_df.iterrows():
        texte = f"{row['Procédé']}\n(Cat: {row['Catalyseur']})"
        plt.annotate(texte, 
                     (row["MSP (€/kg)"], row["GWP_kgCO2_per_kgCPL"]),
                     textcoords="offset points", 
                     xytext=(10, 5), 
                     ha='left',
                     fontsize=8,
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cbd5e1", alpha=0.9))

    plt.title("Front de Pareto : Économie (MSP) vs Écologie (GWP)\nLe compromis idéal de la Superstructure", fontsize=14, fontweight='bold')
    plt.xlabel("MSP - Prix de Vente Minimum (€/kg de Caprolactame)", fontsize=12)
    plt.ylabel("GWP - Empreinte Carbone (kg CO2 eq / kg)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()

def lancer_superstructure_globale():
    df, features, scaler, structure = charger_et_preparer_donnees()
    if df is None: return

    print("\n" + "="*80)
    print(" SUPERSTRUCTURE  ")
    print("="*80)

    try:
        purete_pa6 = float(input("\n -> What is the purity of the incoming plastic waste (%) ? [10-100] : "))
    except ValueError:
        purete_pa6 = 60.0

    try:
        iter_per_proc = int(input(" -> How many initial exploration scenarios per process (ex: 1000) ? : "))
    except ValueError:
        iter_per_proc = 1000

    procs = list(structure.keys())
    tous_les_gagnants_mondiaux = []

    for proc in procs:
        print("\n" + "x"*70)
        print(f"  DEPLOYMENT OF THE DIGITAL TWIN: {proc.replace('_', ' ').upper()}")
        print("x"*70)
        
        expert = entrainer_expert(proc, df, features)
        if expert is None: continue
        
        catalyseurs_possibles = structure[proc]
        scenarios_explores = []

        # --- PHASE 1 : EXPLORATION ---
        print(f"\n   [1/3] Random Exploration Under Physical Constraints ({iter_per_proc} scenarios)...")
        for i in range(iter_per_proc):
            if i > 0 and i % 100 == 0:
                sys.stdout.write(f"\r      Progress : {i} / {iter_per_proc}")
                sys.stdout.flush()

            cat_choisi = str(random.choice(catalyseurs_possibles))
            reacteur_choisi = str(random.choice(REACTEURS_BIBLIO.get(cat_choisi, ["Batch"])))
            
            if "alcoholysis" in proc.lower() or "glycolysis" in proc.lower():
                solvant_choisi = str(random.choice(ALCOOLS_AUTOCLAVE)) if "Autoclave" in reacteur_choisi else str(random.choice(ALCOOLS_BATCH))
            else:
                solvant_choisi = "Water"

            t_min, t_max, p_min_base, p_max_base, is_liquid = obtenir_limites_physiques(proc, cat_choisi, reacteur_choisi)
            temp = random.uniform(t_min, t_max)
            
            if is_liquid:
                p_sat_bar = math.exp(11.68 - 3816.44 / (temp + 227.02))
                p_min = (p_sat_bar * 100) + 200  
                p_max = (p_sat_bar * 100) + 5000 
            else:
                p_min, p_max = p_min_base, p_max_base
                
            press = random.uniform(p_min, p_max)
            time = random.uniform(0.5, 6.0)
            cat_val = 0.0 if cat_choisi == "0" else random.uniform(0.1, 10.0)
            liq_ratio = random.uniform(0.5, 5.0)

            ligne = creer_ligne_ia(features, proc, cat_choisi, reacteur_choisi, temp, press, time, cat_val, solvant_choisi, "0", 0.0, liq_ratio)
            df_test = pd.DataFrame([ligne], columns=features)
            scaled = scaler.transform(df_test)
            rendement = expert.predict(pd.DataFrame(scaled, columns=features))[0]
            
            scenarios_explores.append({
                "Procédé": proc.replace('_', ' '), "Catalyseur": cat_choisi, "Réacteur": reacteur_choisi,
                "Solvant": solvant_choisi, "Temp (°C)": temp, "Temps (h)": time, "Cat (wt%)": cat_val,
                "Pression (KPa)": press, "Ratio Liq": liq_ratio, "Rendement IA (%)": rendement,
                "t_min": t_min, "t_max": t_max, "p_min_base": p_min_base, "p_max_base": p_max_base, "is_liquid": is_liquid
            })
        print(f"\r      Progression : {iter_per_proc} / {iter_per_proc} [Terminé]")

        # --- PHASE 2 : MUTATION GÉNÉTIQUE ---
        print(f"   [2/3] Local Genetic Optimization...")
        top_10_local = sorted(scenarios_explores, key=lambda x: x["Rendement IA (%)"], reverse=True)[:10]
        mutants_locaux = []
        
        for parent in top_10_local:
            for m in range(50):
                new_temp = max(parent["t_min"], min(parent["t_max"], parent["Temp (°C)"] + random.uniform(-10.0, 10.0)))
                
                if parent["is_liquid"]:
                    p_sat_bar = math.exp(11.68 - 3816.44 / (new_temp + 227.02))
                    p_min = (p_sat_bar * 100) + 200
                    p_max = (p_sat_bar * 100) + 5000
                else:
                    p_min, p_max = parent["p_min_base"], parent["p_max_base"]
                    
                new_press = max(p_min, min(p_max, parent["Pression (KPa)"] + random.uniform(-500.0, 500.0)))
                new_time = max(0.1, min(6.0, parent["Temps (h)"] + random.uniform(-0.5, 0.5)))
                new_cat_val = 0.0 if parent["Catalyseur"] == "0" else max(0.1, min(10.0, parent["Cat (wt%)"] + random.uniform(-1.0, 1.0)))
                new_liq_ratio = max(0.5, min(5.0, parent["Ratio Liq"] + random.uniform(-0.5, 0.5)))
                
                ligne = creer_ligne_ia(features, proc, parent["Catalyseur"], parent["Réacteur"], new_temp, new_press, new_time, new_cat_val, parent["Solvant"], "0", 0.0, new_liq_ratio)
                df_test = pd.DataFrame([ligne], columns=features)
                scaled = scaler.transform(df_test)
                rendement_mutant = expert.predict(pd.DataFrame(scaled, columns=features))[0]
                
                mutants_locaux.append({
                    "Procédé": parent["Procédé"], "Catalyseur": parent["Catalyseur"], "Réacteur": parent["Réacteur"],
                    "Solvant": parent["Solvant"], "Temp (°C)": new_temp, "Temps (h)": new_time, "Cat (wt%)": new_cat_val,
                    "Pression (KPa)": new_press, "Ratio Liq": new_liq_ratio, "Rendement IA (%)": rendement_mutant
                })
                
        tous_les_candidats = top_10_local + mutants_locaux
        for c in tous_les_candidats:
            c.pop("t_min", None)
            c.pop("t_max", None)
            c.pop("p_min_base", None)
            c.pop("p_max_base", None)
            c.pop("is_liquid", None)
            
        elite_locale = sorted(tous_les_candidats, key=lambda x: x["Rendement IA (%)"], reverse=True)[:15]

        # --- PHASE 3 : ÉVALUATION BIOSTEAM ---
        print(f"   [3/3] Technical and Economic Evaluation of 15 finalists (Kinetic Mode)...")
        resultats_biosteam_locaux = []
        
        for idx_scenario, scenario in enumerate(elite_locale):
            sys.stdout.write(f"\r      Calculs d'usine : {idx_scenario + 1} / 15 ")
            sys.stdout.flush()

            masse_dechets_entrants = 1000.0 
            masse_nylon_pur = masse_dechets_entrants * (purete_pa6 / 100.0)
            masse_fibres = masse_dechets_entrants - masse_nylon_pur 

            masse_liquide = masse_nylon_pur * scenario["Ratio Liq"]
            masse_cat = masse_nylon_pur * (scenario["Cat (wt%)"] / 100.0)
            
            flux_entrant_brut = {
                "PA6_Solid": masse_nylon_pur,
                "Fibre_Verre": masse_fibres,
                "Caprolactam": 0.0,
                "Oligomers": 0.0
            }

            if scenario["Solvant"] != "Water":
                flux_entrant_brut["Isopropanol"] = masse_liquide
            else:
                flux_entrant_brut["Water"] = masse_liquide

            cat_name = scenario["Catalyseur"]
            if cat_name in ["H3PO4", "HCl", "H2SO4", "MSA", "p-TSA", "BSA", "H3PW12O40 HPA"]:
                flux_entrant_brut["H3PO4_Acide"] = masse_cat
            elif cat_name in ["NaOH", "NaOH-KOH"]:
                flux_entrant_brut["NaOH_Base"] = masse_cat
            elif cat_name in ["m-ZrO2", "0"]:
                flux_entrant_brut["m_ZrO2"] = masse_cat
            elif cat_name in ["Fe- MOF-1", "Fe- MOF-2"]:
                flux_entrant_brut["Fe_MOF_Solid"] = masse_cat

            try:
                kpi = orchestrateur_usine(
                    flux_entrant_brut, scenario["Catalyseur"], scenario["Réacteur"], 
                    scenario["Temp (°C)"], scenario["Pression (KPa)"], 
                    scenario["Temps (h)"], scenario["Solvant"]
                )
                
                production_cpl = kpi.get('Production_CPL_kg_hr', 0)
                
                if production_cpl <= 0:
                    continue
                    
                revenu = production_cpl * 2.0 
                
                cout_vapeur = kpi.get('Cout_Vapeur_hr',0)
                cout_eau = kpi.get('Cout_Eau_Froide_hr',0)
                cout_solvant = kpi.get('Cout_Solvant_hr',0)
                cout_equipement = kpi.get('Cout_Equipement_hr',0)
                cout_reactifs = kpi.get('Cout_Reactifs_hr',0)
                cout_catalyseur = kpi.get('Cout_Catalyseur_hr',0)

                if kpi.get('Type_Usine', "") in ["Neutre", "Neutre_Cinetique"]:
                    cout_dechets = kpi.get('Dechets_Totaux_kg', 0) * 0.002 
                else:
                    cout_dechets = (kpi.get('Eaux_Usees_kg', 0) * 0.002) + (kpi.get('Dechets_Toxiques_kg', 0) * 0.50) 

                FACTEUR_MARGE = 1.45 
                couts_totaux = cout_vapeur + cout_eau + cout_solvant + cout_equipement + cout_reactifs + cout_catalyseur + cout_dechets
            
                msp = (couts_totaux * FACTEUR_MARGE) / production_cpl
                
                scenario["Profit (€/h)"] = revenu - couts_totaux
                scenario["MSP (€/kg)"] = msp
                scenario["Rendement Physique (%)"] = (production_cpl / masse_nylon_pur) * 100.0
                scenario["GWP_kgCO2_per_kgCPL"] = kpi.get('GWP_kgCO2_per_kgCPL', 0.0)
                scenario["Volume (m3)"] = kpi.get("Volume_Reacteur_m3", 0.0)
                scenario["Purete Entrante (%)"] = purete_pa6
                
                resultats_biosteam_locaux.append(scenario)
                
            except Exception as e:
                pass

        print(f"\r      Factory Calculations: 15 / 15 [Completed]       ")   
        
        top_5_local_profit = sorted(resultats_biosteam_locaux, key=lambda x: x.get("Profit (€/h)", -9999), reverse=True)[:5]
        tous_les_gagnants_mondiaux.extend(top_5_local_profit)

    # =====================================================================
    # PHASE 4 : MENU INTERACTIF (BATAILLE ROYALE MONDIALE & PARETO)
    # =====================================================================
    while True:
        print("\n" + "="*90)
        print("  STRATÉGIE D'OPTIMISATION (MENU INTERACTIF)")
        print(" [1] Maximizing Profit / Minimizing MSP (Economics)")
        print(" [2] Minimizing the Carbon Footprint (Environmental)")
        print(" [3] Afficher le Diagramme de Pareto (MSP vs GWP)")
        print(" [4] Quitter le programme")
        choix_opti = input(" -> Your Choice (1/2/3/4) : ")
        
        if choix_opti == "4":
            print("\n Fermeture de la Superstructure. Beau travail !")
            break
            
        elif choix_opti == "3":
            tracer_front_pareto(tous_les_gagnants_mondiaux)
            continue
            
        elif choix_opti == "2":
            classement_final = sorted(tous_les_gagnants_mondiaux, key=lambda x: x.get("GWP_kgCO2_per_kgCPL", 999))
            titre = "ENVIRONMENTAL STRATEGY (MINIMUM GWP)"
            
        elif choix_opti == "1":
            classement_final = sorted(tous_les_gagnants_mondiaux, key=lambda x: x.get("Profit (€/h)", -9999), reverse=True)
            titre = "ECONOMIC STRATEGY (MAXIMUM PROFIT)"
            
        else:
            print(" ❌ Choix invalide. Veuillez taper 1, 2, 3 ou 4.")
            continue

        print(f"\n {titre}")
        print("="*90)
        
        top_5_mondial = classement_final[:5]
        
        for i, res in enumerate(top_5_mondial):
            print(f"\n RANG {i+1} : ** {res['Procédé'].upper()} **")
            print(f"   ->  MSP (Minimum Selling Price) : {res['MSP (€/kg)']:.2f} €/kg de CPL produit")
            print(f"   ->  market profit: {res['Profit (€/h)']:.2f} €/h")
            print(f"   ->  Carbon footprint (GWP)     : {res['GWP_kgCO2_per_kgCPL']:.2f} kg CO2-eq / kg CPL")
            print(f"   ->  Yield (IA predicted)       : {res['Rendement IA (%)']:.2f} %")
            print(f"   ->  Yield (Real Kinetic)       : {res.get('Rendement Physique (%)', 0.0):.2f} %")
            print(f"   ->  Equipment                  : Réacteur = {res['Réacteur']} (Volume calculé : {res.get('Volume (m3)', 0.0):.1f} m³)")
            print(f"   ->  chemistry                  : Catalyseur = {res['Catalyseur']} | Solvant = {res['Solvant']}")
            print(f"   ->  experimental conditions    : {res['Temp (°C)']:.1f} °C | {res['Pression (KPa)']:.0f} KPa | {res['Temps (h)']:.1f} h | Cat: {res['Cat (wt%)']:.1f} % | Ratio Liq: {res['Ratio Liq']:.1f}")
            print("-" * 80)
                
        fichier_sauvegarde = os.path.join(DOSSIER_TRAVAIL, 'superstructure_resultats.csv')
        df_sauvegarde = pd.DataFrame(top_5_mondial)
        df_sauvegarde = df_sauvegarde.round(2) 
        df_sauvegarde.to_csv(fichier_sauvegarde, mode='a', header=not os.path.exists(fichier_sauvegarde), index=False, sep=';', encoding='utf-8')
        print(f"\n  SUCCESS : The ranking has been saved in: {fichier_sauvegarde}")

        # =====================================================================
        # PHASE 5 : TRACÉ DE LA COURBE RENDEMENT VS TEMPS (IA vs CINÉTIQUE)
        # =====================================================================
        if top_5_mondial:
            afficher_courbe = input("\n -> Voulez-vous générer la courbe cinétique du TOP 1 ? (y/n) : ")
            if afficher_courbe.lower() == 'y':
                top1 = top_5_mondial[0]
                print(f"\n--- Génération de la courbe cinétique pour le TOP 1 : {top1['Procédé']} ---")
                
                proc = top1["Procédé"].replace(' ', '_')
                cat = top1["Catalyseur"]
                reacteur = top1["Réacteur"]
                solvant = top1["Solvant"]
                temp = top1["Temp (°C)"]
                press = top1["Pression (KPa)"]
                cat_wt = top1["Cat (wt%)"]
                liq_ratio = top1["Ratio Liq"]
                
                if "Acidic" in proc: proc_type = "acidic"
                elif "Neutral" in proc: proc_type = "neutral"
                elif "Alkaline" in proc: proc_type = "alkaline"
                else: proc_type = "alcoholysis"
                
                params_cinetiques = obtenir_params_pour_reacteur(proc_type, cat, solvant)
                Ea = params_cinetiques.get("Ea_J_mol", 105000.0)
                A = params_cinetiques.get("A_h_1", 3.5e9)
                n = params_cinetiques.get("n_ordre", 1.0)
                m = params_cinetiques.get("m_catalyseur", 0.0)
                
                temps_heures = np.linspace(0.1, 10.0, 50)
                rendements_ia = []
                rendements_cinetiques = []
                
                expert = entrainer_expert(proc, df, features)
                
                for t in temps_heures:
                    ligne = creer_ligne_ia(features, proc, cat, reacteur, temp, press, t, cat_wt, solvant, "0", 0.0, liq_ratio)
                    df_test = pd.DataFrame([ligne], columns=features)
                    scaled = scaler.transform(df_test)
                    r_ia = expert.predict(pd.DataFrame(scaled, columns=features))[0]
                    rendements_ia.append(np.clip(r_ia, 0.0, 100.0))
                    
                    T_k = temp + 273.15
                    k = A * np.exp(-Ea / (8.314 * T_k)) * (max(cat_wt, 1e-5)**m)
                    
                    if abs(n - 1.0) < 1e-4:
                        X = 1.0 - np.exp(-k * t)
                    else:
                        base = 1.0 + k * (n - 1.0) * t
                        X = 1.0 - (max(base, 1e-10) ** (1.0 / (1.0 - n)))
                    
                    X = np.clip(X, 0.0, 0.999)
                    rendements_cinetiques.append(X * 100.0)
                    
                plt.figure(figsize=(10, 6))
                plt.plot(temps_heures, rendements_ia, label="Modèle IA (Random Forest)", color="#1f77b4", linewidth=2.5)
                plt.plot(temps_heures, rendements_cinetiques, label="Modèle Physique (Cinétique Arrhenius)", color="#d62728", linestyle="--", linewidth=2.5)
                
                t_opt = top1["Temps (h)"]
                y_opt_ia = top1["Rendement IA (%)"]
                y_opt_kin = top1.get("Rendement Physique (%)", 0.0)
                
                plt.scatter([t_opt], [y_opt_ia], color="#1f77b4", s=100, zorder=5, edgecolors='black')
                plt.scatter([t_opt], [y_opt_kin], color="#d62728", s=100, zorder=5, edgecolors='black')
                plt.axvline(x=t_opt, color='gray', linestyle=':', label=f"Temps de séjour choisi ({t_opt:.1f} h)")
                
                plt.title(f"Dynamique Temporelle du TOP 1 : {top1['Procédé']}\nCatalyseur: {cat} | T = {temp:.1f}°C | P = {press:.0f} kPa | Cat = {cat_wt:.1f} wt%", fontsize=12, pad=15)
                plt.xlabel("Temps de réaction (heures)", fontsize=11)
                plt.ylabel("Rendement en Caprolactame (%)", fontsize=11)
                plt.xlim(0, 10)
                plt.ylim(0, 105)
                plt.legend(loc="lower right", fontsize=10)
                plt.grid(True, linestyle="--", alpha=0.6)
                plt.tight_layout()
                
                chemin_graphe = os.path.join(DOSSIER_TRAVAIL, 'top1_cinetique_vs_ia.png')
                plt.savefig(chemin_graphe, dpi=300)
                print(f"  -> Graphique sauvegardé avec succès dans : {chemin_graphe}")
                
                plt.show()

if __name__ == "__main__":
    lancer_superstructure_globale()