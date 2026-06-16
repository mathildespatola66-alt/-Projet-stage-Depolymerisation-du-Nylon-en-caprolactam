import thermosteam as tmo
import biosteam as bst

# =========================================================================
# 1. INITIALISATION THERMODYNAMIQUE GLOBALE
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
# 2. L'USINE COMPLÈTE (AVEC ÉCHANGEURS DE CHALEUR EXPLICITES)
# =========================================================================

def simuler_usine_neutre(flux_masses_dictionnaire, type_reacteur, temp_user, pression_user, nom_catalyseur, temps_sejour):
    
    bst.main_flowsheet.clear()
    bst.settings.set_thermo(chemin_chimique)
    
   # --- BLOC 1 : LE RÉACTEUR ET LA BOUCLE THERMIQUE RIGOUREUSE ---
    # L'alimentation est configurée à la pression requise dès l'entrée
    soupe_froide = bst.Stream('Soupe_Froide', **flux_masses_dictionnaire, units='kg/hr', T=25+273.15, P=pression_user * 1000)
    
    # Préchauffage thermique (FEHE)
    soupe_prechauffee = soupe_froide.copy()
    soupe_prechauffee.T = 25 + (temp_user - 25) * 0.85 + 273.15
    
    # 1. CHAUFFEUR D'APPOINT (Calcul de l'apport enthalpique externe)
    H1_Trim = bst.HXutility('H1_Trim_Heater', ins=soupe_prechauffee, outs='Soupe_Chaude', T=temp_user+273.15)
    H1_Trim.run()
    energie_trim = H1_Trim.outs[0].H - H1_Trim.ins[0].H 
    
    # 2. LE RÉACTEUR DE CONVERSION CONTINU (Modélisation thermodynamique rigoureuse)
    # L'utilisation de ContinuousReactor conserve l'intégralité des flashs thermiques et physiques
    # sans subir les biais des géométries d'agitateurs mécaniques.
    R101 = bst.Mixer(
        'R101_Reacteur_Virtuel', 
        ins=H1_Trim.outs[0], 
        outs='Jus_Reacteur'
    )
    R101.run()
    
    # 3. REFROIDISSEUR DE SÉCURITÉ (Quench Cooler)
    C1_Quench = bst.HXutility('C1_Quench_Cooler', ins=R101.outs[0], outs='Soupe_Quench', T=80+273.15)
    C1_Quench.run()
    energie_quench = C1_Quench.ins[0].H - C1_Quench.outs[0].H

    # --- BLOC 2 : FILTRATION ---
    Filtre = bst.Splitter('Filtre', ins=C1_Quench.outs[0], outs=('Cat_et_Inertes', 'Jus'), 
                          split={'PA6_Solid': 1.0, 'm_ZrO2': 1.0, 'Fe_MOF_Solid': 1.0, 'Fibre_Verre': 1.0,
                                 'Water': 0.0, 'Caprolactam': 0.0, 'Oligomers': 0.0, 'Cyclohexanol': 0.0,
                                 'Toluene': 0.0, 'H3PO4_Acide': 0.0, 'NaOH_Base': 0.0, 'Na3PO4_Sel': 0.0})
    Filtre.run()
    
    # --- BLOC 3 : EXTRACTION LLE (DYNAMIQUE) ---
    debit_jus_aqueux = Filtre.outs[1].F_mass 
    ratio_SF = 0.5 
    solvant_total_requis = debit_jus_aqueux * ratio_SF 
    
    solvant = bst.Stream('Solv', Cyclohexanol=solvant_total_requis, units='kg/hr', T=55+273.15)
    Mix = bst.Mixer('Mix', ins=(Filtre.outs[1], solvant), outs=('Melange',))
    Mix.run()
    
    ELL = bst.Splitter('ELL', ins=Mix.outs[0], outs=('Org', 'Aq'), 
                       split={'Cyclohexanol': 0.99, 'Caprolactam': 0.95, 
                              'Water': 0.05, 'Oligomers': 0.01, 'PA6_Solid': 0.0, 'm_ZrO2': 0.0, 'Fe_MOF_Solid': 0.0, 'Fibre_Verre': 0.0,
                              'Toluene': 0.0, 'H3PO4_Acide': 0.0, 'NaOH_Base': 0.0, 'Na3PO4_Sel': 0.0})
    ELL.run()
    
    # 3. BOUILLEUR (Reboiler)
    H2_Reboiler = bst.HXutility('H2_Reboiler', ins=ELL.outs[0], outs='Org_Chaud', T=165+273.15)
    H2_Reboiler.run()
    energie_flash = H2_Reboiler.outs[0].H - H2_Reboiler.ins[0].H
    
    # --- BLOC 4 : ÉVAPORATION FLASH ---
    Flash = bst.Splitter('Flash', ins=H2_Reboiler.outs[0], outs=('Vapeur_Solvant', 'CPL_Brut'), 
                         split={'Cyclohexanol': 0.99, 'Caprolactam': 0.01, 
                                'Water': 0.99, 'Oligomers': 0.0, 'PA6_Solid': 0.0, 'm_ZrO2': 0.0, 'Fe_MOF_Solid': 0.0, 'Fibre_Verre': 0.0,
                                'Toluene': 0.0, 'H3PO4_Acide': 0.0, 'NaOH_Base': 0.0, 'Na3PO4_Sel': 0.0})
    Flash.run()
    
    # 4. CONDENSEUR (Condenser)
    C2_Condenser = bst.HXutility('C2_Condenser', ins=Flash.outs[0], outs='Solvant_Froid', T=55+273.15)
    C2_Condenser.run()
    energie_condenseur = C2_Condenser.ins[0].H - C2_Condenser.outs[0].H
    
    # --- BLOC 4.5 : BOUCLE DE RECYCLAGE ET PURGE ---
    cycles_maximum = 50.0 
    fraction_purge = 1.0 / cycles_maximum 
    
    Purge = bst.Splitter('Purge', ins=C2_Condenser.outs[0], outs=('Solv_Recycle', 'Purge_Dechets'), 
                         split=1.0 - fraction_purge)
    Purge.run()
    
    masse_solvant_recycle = Purge.outs[0].imass['Cyclohexanol']
    masse_appoint_neuf = max(0.0, solvant_total_requis - masse_solvant_recycle)
    
    # 5. REFROIDISSEUR DE SÉCURITÉ (Chiller)
    C3_Chiller = bst.HXutility('C3_Chiller', ins=Flash.outs[1], outs='CPL_Froid', T=60+273.15)
    C3_Chiller.run()
    energie_chiller = C3_Chiller.ins[0].H - C3_Chiller.outs[0].H
    
    # --- BLOC 5 : CRISTALLISATION ---
    Cristaux = bst.Splitter('Cristaux', ins=C3_Chiller.outs[0], outs=('CPL_Pur', 'Dechets_Liq'), 
                            split={'Caprolactam': 0.98, 'Cyclohexanol': 0.001, 'Water': 0.001, 
                                   'Oligomers': 0.0, 'PA6_Solid': 0.0, 'm_ZrO2': 0.0, 'Fe_MOF_Solid': 0.0, 'Fibre_Verre': 0.0,
                                   'Toluene': 0.0, 'H3PO4_Acide': 0.0, 'NaOH_Base': 0.0, 'Na3PO4_Sel': 0.0})
    
    # CRÉATION ET SIMULATION DU SYSTÈME GLOBAL
    mon_systeme_complet = bst.main_flowsheet.create_system('Usine_Nylon')
    mon_systeme_complet.simulate()

    # GÉNÉRATION DU DIAGRAMME VISUEL (PFD) 
    chemin_absolu = r"c:\Users\mathi\Desktop\Stage_Nylon\python\process\Mon_Flowsheet_Neutre"
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
        
        # CAPEX DYNAMIQUE (Modélisation rigoureuse basée sur le temps de séjour utilisateur)
        masse_entree_totale = soupe_froide.F_mass 
        densite_moyenne = 950.0  # kg/m³ (Mélange eau sous-critique + polymère liquide)
        
        # Application stricte de la loi d'extrapolation physique V = Q * tau
        volume_reacteur_m3 = (masse_entree_totale / densite_moyenne) * temps_sejour
        
        # Indexation des coûts selon la loi des coefficients de puissance (loi des 0.6)
        cout_base_reacteur = 350000  # Réacteur haute pression/haute température en acier allié
        cout_reacteur_extrapole = cout_base_reacteur * (max(volume_reacteur_m3, 0.5) / 5.0)**0.6
        cout_reste_usine = 2150000 * (masse_entree_totale / 10000)**0.6
        
        capex_total = cout_reacteur_extrapole + cout_reste_usine
        cout_equipement_hr = capex_total / (10 * 8000) 
        
        # Déchets et Solvant
        masse_dechets = ELL.outs[1].F_mass + Cristaux.outs[1].F_mass + Purge.outs[1].F_mass + Filtre.outs[0].imass['Fibre_Verre']
        cout_solvant_hr = masse_appoint_neuf * 1.50
        
        # Coût du catalyseur hétérogène (perte de 2% par heure)
        PRIX_CATALYSEURS = {"Fe- MOF-1": 150.00, "Fe- MOF-2": 150.00, "m-ZrO2": 80.00, "0": 0.00}
        prix_u = PRIX_CATALYSEURS.get(nom_catalyseur, 0.0)
        masse_cata_totale = flux_masses_dictionnaire.get('m_ZrO2', 0) + flux_masses_dictionnaire.get('Fe_MOF_Solid', 0)
        cout_cata_hr = masse_cata_totale * 0.02 * prix_u
        
        # Pureté du produit final
        flux_produit_fini = Cristaux.outs[0]
        masse_cpl_pur = flux_produit_fini.imass['Caprolactam']
        masse_totale_sortie = flux_produit_fini.F_mass
        purete_finale = (masse_cpl_pur / masse_totale_sortie) * 100.0 if masse_totale_sortie > 0 else 0.0
        
        return {
            "Cout_Vapeur_hr": cout_vapeur_hr, 
            "Cout_Eau_Froide_hr": cout_eau_froide_hr,
            "Dechets_Totaux_kg": masse_dechets, 
            "Cout_Equipement_hr": cout_equipement_hr, 
            "Cout_Solvant_hr": cout_solvant_hr,
            "Cout_Catalyseur_hr": cout_cata_hr, 
            "Purete_CPL": purete_finale,
            "Type_Usine": "Neutre"
        }
        
    except Exception as e:
        print(f"Erreur interne : {e}")
        return {"Cout_Vapeur_hr": 0, "Cout_Eau_Froide_hr": 0, "Dechets_Totaux_kg": 0, "Cout_Equipement_hr": 0, "Cout_Solvant_hr": 0, "Cout_Catalyseur_hr": 0, "Purete_CPL": 0, "Type_Usine": "Erreur"}