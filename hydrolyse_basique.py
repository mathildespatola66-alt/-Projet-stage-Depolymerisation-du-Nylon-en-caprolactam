import thermosteam as tmo
import biosteam as bst

# =========================================================================
# 1. INITIALISATION THERMODYNAMIQUE GLOBALE (À METTRE PARTOUT)
# =========================================================================
Water = bst.Chemical('Water')
CPL = bst.Chemical('Caprolactam')
Cyclohexanol = bst.Chemical('Cyclohexanol')
Toluene = bst.Chemical('Toluene')

# Création des "clones" physiques pour les solides et les catalyseurs
PA6 = bst.Chemical('Glycerol').copy('PA6_Solid')
Catalyseur_Solide = bst.Chemical('Glycerol').copy('m_ZrO2')
Fe_MOF_Solid = bst.Chemical('Glycerol').copy('Fe_MOF_Solid') 
Oligomers = bst.Chemical('Glycerol').copy('Oligomers')
Fibre_Verre = bst.Chemical('Glycerol').copy('Fibre_Verre')

# Déclaration des acides, bases et sels
H3PO4_Acide = bst.Chemical('Glycerol').copy('H3PO4_Acide')
NaOH_Base = bst.Chemical('Glycerol').copy('NaOH_Base')
Na3PO4_Sel = bst.Chemical('Glycerol').copy('Na3PO4_Sel')

# On rassemble TOUT LE MONDE dans le dictionnaire
chemin_chimique = bst.Chemicals([
    Water, CPL, Cyclohexanol, Toluene, Oligomers, PA6, 
    Catalyseur_Solide, Fe_MOF_Solid, Fibre_Verre, 
    H3PO4_Acide, NaOH_Base, Na3PO4_Sel
])
bst.settings.set_thermo(chemin_chimique)
CPL.smiles = 'C1CCC(=O)NCC1'

# =========================================================================
# 2. L'USINE COMPLÈTE (MOTEUR THERMODYNAMIQUE)
# =========================================================================
def simuler_usine_basique(flux_masses_dictionnaire, type_reacteur, temp_user, pression_user, nom_catalyseur, temps_sejour):
    bst.main_flowsheet.clear() 
    bst.settings.set_thermo(chemin_chimique)
    
    # --- BLOC 1 : LE RÉACTEUR ET LA BOUCLE THERMIQUE ---
    # CORRECTION : Entrée sous haute pression
    soupe_froide = bst.Stream('Soupe_Froide', **flux_masses_dictionnaire, units='kg/hr', T=25+273.15, P=pression_user * 1000)
    
    soupe_prechauffee = soupe_froide.copy()
    soupe_prechauffee.T = 25 + (temp_user - 25) * 0.85 + 273.15
    
    H1_Trim = bst.HXutility('H1_Trim_Heater', ins=soupe_prechauffee, outs='Soupe_Chaude', T=temp_user+273.15)
    H1_Trim.run()
    energie_trim = H1_Trim.outs[0].H - H1_Trim.ins[0].H 
    
    # CORRECTION : Réacteur macroscopique rigoureux
    R101 = bst.Mixer(
        'R101_Reacteur_Basique', 
        ins=H1_Trim.outs[0], 
        outs='Jus_Reacteur'
    )
    R101.run()
    
    # CORRECTION : Entrée corrigée pour le Quench
    C1_Quench = bst.HXutility('C1_Quench_Cooler', ins=R101.outs[0], outs='Soupe_Quench', T=80+273.15)
    C1_Quench.run()
    energie_quench = C1_Quench.ins[0].H - C1_Quench.outs[0].H 
    
    # --- BLOC 1.5 : NEUTRALISATION (Spécifique voie Basique) ---
    masse_base_presente = flux_masses_dictionnaire.get('NaOH_Base', 0)
    # On achète de l'acide (H3PO4) pour détruire la soude (ratio de sécurité 1.5)
    masse_acide_achetee = masse_base_presente * 1.5
    flux_acide = bst.Stream('Flux_Acide_Achete', H3PO4_Acide=masse_acide_achetee, units='kg/hr', T=25+273.15)
    
    Cuve_Neutralisation = bst.Mixer('Cuve_Neutralisation', ins=(C1_Quench.outs[0], flux_acide), outs=('Jus_Neutralise',))
    Cuve_Neutralisation.run()
    
    # --- BLOC 2 : FILTRATION DES SELS ET INERTES ---
    Filtre_Sels = bst.Splitter('Filtre_Sels', ins=Cuve_Neutralisation.outs[0], outs=('Dechets_Solides', 'Jus_Clair'), 
                               split={'PA6_Solid': 1.0, 'Fibre_Verre': 1.0, 'NaOH_Base': 1.0, 'H3PO4_Acide': 1.0, 'Na3PO4_Sel': 1.0,
                                      'Water': 0.0, 'Caprolactam': 0.0, 'Oligomers': 0.0, 'Toluene': 0.0})
    Filtre_Sels.run()
    
    # --- BLOC 3 : EXTRACTION LLE (DYNAMIQUE AU TOLUÈNE) ---
    debit_jus_aqueux = Filtre_Sels.outs[1].F_mass 
    ratio_SF = 0.5 
    solvant_total_requis = debit_jus_aqueux * ratio_SF 
    
    solvant = bst.Stream('Solv_Toluene', Toluene=solvant_total_requis, units='kg/hr', T=55+273.15)
    Mix = bst.Mixer('Mix', ins=(Filtre_Sels.outs[1], solvant), outs=('Melange',))
    Mix.run()
    
    ELL = bst.Splitter('ELL', ins=Mix.outs[0], outs=('Org', 'Aq'), 
                       split={'Toluene': 0.99, 'Caprolactam': 0.95, 
                              'Water': 0.05, 'Oligomers': 0.01, 'PA6_Solid': 0.0, 'Fibre_Verre': 0.0, 'NaOH_Base': 0.0, 'H3PO4_Acide': 0.0})
    ELL.run()
    
    # --- BLOC 4 : ÉVAPORATION FLASH ---
    H2_Reboiler = bst.HXutility('H2_Reboiler', ins=ELL.outs[0], outs='Org_Chaud', T=150+273.15)
    H2_Reboiler.run()
    energie_flash = H2_Reboiler.outs[0].H - H2_Reboiler.ins[0].H
    
    Flash = bst.Splitter('Flash', ins=H2_Reboiler.outs[0], outs=('Vapeur_Solvant', 'CPL_Brut'), 
                         split={'Toluene': 0.99, 'Caprolactam': 0.01, 
                                'Water': 0.99, 'Oligomers': 0.0, 'PA6_Solid': 0.0, 'Fibre_Verre': 0.0})
    Flash.run()
    
    # CONDENSEUR DE TOLUÈNE
    C2_Condenser = bst.HXutility('C2_Condenser', ins=Flash.outs[0], outs='Solvant_Froid', T=55+273.15)
    C2_Condenser.run()
    energie_condenseur = C2_Condenser.ins[0].H - C2_Condenser.outs[0].H
    
    # --- BLOC 4.5 : BOUCLE DE RECYCLAGE ET PURGE ---
    cycles_maximum = 50.0 
    fraction_purge = 1.0 / cycles_maximum 
    
    Purge = bst.Splitter('Purge', ins=C2_Condenser.outs[0], outs=('Solv_Recycle', 'Purge_Dechets'), split=1.0 - fraction_purge)
    Purge.run()
    
    masse_solvant_recycle = Purge.outs[0].imass['Toluene']
    masse_appoint_neuf = max(0.0, solvant_total_requis - masse_solvant_recycle)
    
    # --- BLOC 5 : CRISTALLISATION ---
    C3_Chiller = bst.HXutility('C3_Chiller', ins=Flash.outs[1], outs='CPL_Froid', T=60+273.15)
    C3_Chiller.run()
    energie_chiller = C3_Chiller.ins[0].H - C3_Chiller.outs[0].H
    
    Cristaux = bst.Splitter('Cristaux', ins=C3_Chiller.outs[0], outs=('CPL_Pur', 'Dechets_Liq'), 
                            split={'Caprolactam': 0.98, 'Toluene': 0.001, 'Water': 0.001, 
                                   'Oligomers': 0.0, 'PA6_Solid': 0.0, 'm_ZrO2': 0.0, 'Fe_MOF_Solid': 0.0, 'Fibre_Verre': 0.0})
    
    # CRÉATION DU SYSTÈME GLOBAL
    mon_systeme_complet = bst.main_flowsheet.create_system('Usine_Basique')
    mon_systeme_complet.simulate()

    # GÉNÉRATION DU NOUVEAU DIAGRAMME VISUEL (PFD) 
    chemin_absolu = r"c:\Users\mathi\Desktop\Stage_Nylon\python\process\Mon_Flowsheet_Basique"
    mon_systeme_complet.diagram(format='png', file=chemin_absolu)

    # =====================================================================
    # 3. CALCUL DES INDICATEURS DE PERFORMANCE (KPI)
    # =====================================================================
    try:
        # Énergies
        energie_chaude_totale_kJ = abs(energie_trim) + abs(energie_flash)
        energie_froide_totale_kJ = abs(energie_quench) + abs(energie_condenseur) + abs(energie_chiller)
        cout_vapeur_hr = (energie_chaude_totale_kJ / 1000) * 0.005
        cout_eau_froide_hr = (energie_froide_totale_kJ / 1000) * 0.0005
        
        # CORRECTION : CAPEX DYNAMIQUE
        masse_entree_totale = soupe_froide.F_mass 
        densite_moyenne = 950.0  # kg/m³
        volume_reacteur_m3 = (masse_entree_totale / densite_moyenne) * temps_sejour
        
        cout_base_reacteur = 350000 
        cout_reacteur_extrapole = cout_base_reacteur * (max(volume_reacteur_m3, 0.5) / 5.0)**0.6
        cout_reste_usine = 2150000 * (masse_entree_totale / 10000)**0.6
        
        capex_total = cout_reacteur_extrapole + cout_reste_usine
        cout_equipement_hr = capex_total / (10 * 8000) 
        
        # Déchets : Séparation Eau vs Toxique
        masse_eaux_usees = ELL.outs[1].F_mass + Cristaux.outs[1].F_mass + Purge.outs[1].F_mass
        masse_sels_toxiques = Filtre_Sels.outs[0].F_mass
        
        cout_solvant_hr = masse_appoint_neuf * 1.50
        cout_acide_hr = masse_acide_achetee * 0.80 # On achète de l'acide pour neutraliser
        
        # Pureté
        masse_cpl_pur = Cristaux.outs[0].imass['Caprolactam']
        masse_totale_sortie = Cristaux.outs[0].F_mass
        purete_finale = (masse_cpl_pur / masse_totale_sortie) * 100.0 if masse_totale_sortie > 0 else 0.0
        
        PRIX_CATALYSEURS = {"NaOH-KOH": 0.40, "NaOH": 0.30}
        prix_u = PRIX_CATALYSEURS.get(nom_catalyseur, 0.30)
        masse_cata = flux_masses_dictionnaire.get('NaOH_Base', 0)
        
        cout_cata_hr = masse_cata * prix_u # On paie 100% de la soude entrant dans le réacteur
        
        return {
            "Cout_Vapeur_hr": cout_vapeur_hr, 
            "Cout_Eau_Froide_hr": cout_eau_froide_hr,
            "Eaux_Usees_kg": masse_eaux_usees,           
            "Dechets_Toxiques_kg": masse_sels_toxiques,  
            "Cout_Equipement_hr": cout_equipement_hr, 
            "Cout_Solvant_hr": cout_solvant_hr,
            "Cout_Reactifs_hr": cout_acide_hr,
            "Cout_Catalyseur_hr": cout_cata_hr, 
            "Purete_CPL": purete_finale,
            "Type_Usine": "Basique"
        }
        
    except Exception as e:
        print(f"Erreur interne : {e}")
        # Dictionnaire d'erreur aligné avec les vraies clés !
        return {"Cout_Vapeur_hr": 0, "Cout_Eau_Froide_hr": 0, "Eaux_Usees_kg": 0, "Dechets_Toxiques_kg": 0, "Cout_Equipement_hr": 0, "Cout_Solvant_hr": 0, "Cout_Reactifs_hr": 0, "Cout_Catalyseur_hr": 0, "Purete_CPL": 0, "Type_Usine": "Erreur"}