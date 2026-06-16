from hydrolyse_neutre import simuler_usine_neutre
from hydrolyse_acide import simuler_usine_acide
from hydrolyse_basique import simuler_usine_basique

def generer_flux_sortie_reacteur(masse_pa6_in, masse_eau_in, type_catalyseur, masse_catalyseur, predictions_IA, masse_fibres):
    conversion = predictions_IA['Conversion']
    rendement = predictions_IA['Rendement_CPL']
    
    masse_cpl = masse_pa6_in * rendement
    masse_oligomeres = masse_pa6_in * (conversion - rendement)
    masse_pa6_restant = masse_pa6_in * (1 - conversion)
    
    flux = {
        'Water': masse_eau_in,
        'Caprolactam': masse_cpl,
        'Oligomers': masse_oligomeres,
        'PA6_Solid': masse_pa6_restant,
        'Fibre_Verre': masse_fibres
    }
    
    cat_nom = str(type_catalyseur).upper()
    
    # Assignation du catalyseur (On assimile les acides/bases forts aux modèles de base pour la thermo)
    if 'ZRO2' in cat_nom:
        flux['m_ZrO2'] = masse_catalyseur
    elif 'MOF' in cat_nom:
        flux['Fe_MOF_Solid'] = masse_catalyseur
    elif 'H3PO4' in cat_nom or 'HCL' in cat_nom or 'H2SO4' in cat_nom or 'MSA' in cat_nom or 'TSA' in cat_nom or 'BSA' in cat_nom or 'HPA' in cat_nom:
        flux['H3PO4_Acide'] = masse_catalyseur 
    elif 'NAOH' in cat_nom or 'KOH' in cat_nom:
        flux['NaOH_Base'] = masse_catalyseur
        
    return flux

# =========================================================================
# CORRECTION : Ajout de l'argument 'temps_sejour'
# =========================================================================
def orchestrateur_usine(flux_biosteam, type_catalyseur, type_reacteur, temp_user, pression_user, temps_sejour):
    cat_nom = str(type_catalyseur).upper()
    
    if 'ZRO2' in cat_nom or 'MOF' in cat_nom or str(type_catalyseur) == '0':
        # On passe le temps de séjour à l'usine neutre
        return simuler_usine_neutre(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour)
        
    elif 'H3PO4' in cat_nom or 'HCL' in cat_nom or 'H2SO4' in cat_nom or 'MSA' in cat_nom or 'TSA' in cat_nom or 'BSA' in cat_nom or 'HPA' in cat_nom:
        # On passe le temps de séjour à l'usine acide
        return simuler_usine_acide(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour)
        
    elif 'NAOH' in cat_nom or 'KOH' in cat_nom:
        # On passe le temps de séjour à l'usine basique
        return simuler_usine_basique(flux_biosteam, type_reacteur, temp_user, pression_user, type_catalyseur, temps_sejour)
        
    else:
        raise ValueError(f"ERREUR : Catalyseur '{type_catalyseur}' inconnu dans l'orchestrateur !")