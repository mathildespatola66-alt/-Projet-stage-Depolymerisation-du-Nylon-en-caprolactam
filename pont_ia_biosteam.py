import json
import os

from hydrolyse_neutre import simuler_usine_neutre
from hydrolyse_acide import simuler_usine_acide
from hydrolyse_basique import simuler_usine_basique
from alcoolyse import simuler_usine_alcoolyse

def charger_toutes_les_cinetiques():
    chemin = r"c:\Users\mathi\Desktop\Stage_Nylon\python\superstructura\parametres_cinetiques.json"
    try:
        with open(chemin, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def obtenir_params_pour_reacteur(type_process, cat, solvant):
    data = charger_toutes_les_cinetiques()
    cat_nom = str(cat).upper()
    solvant_nom = str(solvant).lower()
    
    fallback = {"Ea_J_mol": 105000.0, "A_h_1": 3.5e9, "n_ordre": 1.0, "m_catalyseur": 0.0}
    
    if "alcoholysis" in type_process.lower():
        if any(x in solvant_nom for x in ["2-propanol", "2-butanol", "cyclohexanol", "cyclododecanol", "2-methyl-2-propanol", "2-methyl-2-butanol"]):
            return data.get("alcoolyse_secondaire", fallback)
        return data.get("alcoolyse_primaire", fallback)
        
    elif "acidic" in type_process.lower():
        # LIAISON DE TOUS LES ACIDES AVEC LE JSON
        if "HPA" in cat_nom: return data.get("acide_hpa", fallback)
        if "MOF-1" in cat_nom: return data.get("acide_mof1", fallback)
        if "MOF-2" in cat_nom: return data.get("acide_mof2", fallback)
        if "H3PO4" in cat_nom: return data.get("acide_h3po4", fallback)
        if "HCl" in cat_nom: return data.get("acide_hcl", fallback)
        return fallback
        
    elif "neutral" in type_process.lower():
        if "m-ZrO2" in cat_nom: return data.get("neutre_m-ZrO2", fallback)
        return data.get("neutre", fallback)
        
    elif "alkaline" in type_process.lower():
        if "NaOH-KOH" in cat_nom: return data.get("alcaline_naoh-koh", fallback)
        return data.get("alcaline", fallback)
        
    return fallback

def orchestrateur_usine(flux_biosteam, type_catalyseur, type_reacteur, temp_user, pression_user, temps_sejour, nom_solvant="Water"):
    cat_nom = str(type_catalyseur).upper()
    
    if nom_solvant != "Water":
        proc_type = "alcoholysis"
        params = obtenir_params_pour_reacteur(proc_type, type_catalyseur, nom_solvant)
        print(f"  [Debug] Procédé: {proc_type} | Catalyseur: {cat_nom} | Ea injectée: {params['Ea_J_mol']} J/mol")
        return simuler_usine_alcoolyse(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour, nom_solvant, params=params)
    
    elif 'ZRO2' in cat_nom or 'MOF' in cat_nom or str(type_catalyseur) == '0':
        proc_type = "neutral"
        params = obtenir_params_pour_reacteur(proc_type, type_catalyseur, nom_solvant)
        print(f"  [Debug] Procédé: {proc_type} | Catalyseur: {cat_nom} | Ea injectée: {params['Ea_J_mol']} J/mol")
        return simuler_usine_neutre(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour, params=params)
        
    # L'ORCHESTRATEUR ROUTE TOUS LES ACIDES VERS L'USINE ACIDE
    elif 'H3PO4' in cat_nom or 'HCL' in cat_nom or 'HPA' in cat_nom:
        proc_type = "acidic"
        params = obtenir_params_pour_reacteur(proc_type, type_catalyseur, nom_solvant)
        print(f"  [Debug] Procédé: {proc_type} | Catalyseur: {cat_nom} | Ea injectée: {params['Ea_J_mol']} J/mol")
        return simuler_usine_acide(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour, params=params)
        
    elif 'NAOH' in cat_nom or 'KOH' in cat_nom:
        proc_type = "alkaline"
        params = obtenir_params_pour_reacteur(proc_type, type_catalyseur, nom_solvant)
        print(f"  [Debug] Procédé: {proc_type} | Catalyseur: {cat_nom} | Ea injectée: {params['Ea_J_mol']} J/mol")
        return simuler_usine_basique(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour, params=params)
        
    else:
        return {"Type_Usine": "Erreur", "Production_CPL_kg_hr": 0.0}
    
def generer_flux_sortie_reacteur(masse_pa6_in, masse_eau_in, type_catalyseur, masse_catalyseur, predictions_IA, masse_fibres, nom_solvant="Water"):
    try:
        conversion = predictions_IA.get('Conversion', 0.5)
        rendement = predictions_IA.get('Rendement_CPL', 0.4)
    except:
        conversion, rendement = 0.5, 0.4
        
    masse_cpl = masse_pa6_in * rendement
    masse_pa6_restant = masse_pa6_in * (1 - conversion)
    
    flux = {'Caprolactam': masse_cpl, 'PA6_Solid': masse_pa6_restant, 'Fibre_Verre': masse_fibres}
    
    if nom_solvant != "Water": 
        flux['Isopropanol'] = masse_eau_in
    else: 
        flux['Water'] = masse_eau_in
        
    return flux