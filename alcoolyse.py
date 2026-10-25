import thermosteam as tmo
import biosteam as bst
import numpy as np

Water = bst.Chemical('Water')
CPL = bst.Chemical('Caprolactam')
Cyclohexanol = bst.Chemical('Cyclohexanol')
Toluene = bst.Chemical('Toluene')
Isopropanol = bst.Chemical('Isopropanol')
PA6 = bst.Chemical('Glycerol').copy('PA6_Solid')
m_ZrO2 = bst.Chemical('Glycerol').copy('m_ZrO2')
Catalyseur_Solide = bst.Chemical('Glycerol').copy('Catalyseur_Solide')
Fe_MOF_Solid = bst.Chemical('Glycerol').copy('Fe_MOF_Solid') 
Oligomers = bst.Chemical('Glycerol').copy('Oligomers')
Fibre_Verre = bst.Chemical('Glycerol').copy('Fibre_Verre')
H3PO4_Acide = bst.Chemical('Glycerol').copy('H3PO4_Acide')
NaOH_Base = bst.Chemical('Glycerol').copy('NaOH_Base')
Na3PO4_Sel = bst.Chemical('Glycerol').copy('Na3PO4_Sel')

chemin_chimique = bst.Chemicals([Water, CPL, Cyclohexanol, Toluene, Isopropanol, Oligomers, PA6, m_ZrO2, Catalyseur_Solide, Fe_MOF_Solid, Fibre_Verre, H3PO4_Acide, NaOH_Base, Na3PO4_Sel])
bst.settings.set_thermo(chemin_chimique)
CPL.smiles = 'C1CCC(=O)NCC1'

class ReacteurCinetiqueGenerique(bst.Unit):
    _ins = ('jus_entrant',)
    _outs = ('jus_sortant',)
    
    def __init__(self, ID='', ins=None, outs=(), T_op=250.0, tau_h=2.0, params=None):
        bst.Unit.__init__(self, ID, ins, outs)
        self.T_op = T_op + 273.15
        self.tau = tau_h
        if params is None: params = {}
        self.Ea = params.get("Ea_J_mol", 100000.0)
        self.A = params.get("A_h_1", 1e9)
        self.n = params.get("n_ordre", 1.0)
        self.m = params.get("m_catalyseur", 0.0)

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.T = self.T_op
        
        cat = 0.0
        masse_tot = max(feed.F_mass, 1e-6)
        
        if feed.imass['Catalyseur_Solide'] > 0: 
            cat = feed.imass['Catalyseur_Solide'] / masse_tot * 100
        
        k = self.A * np.exp(-self.Ea / (8.314 * self.T_op)) * (max(cat, 1e-5)**self.m)
        if abs(self.n - 1.0) < 1e-4: X = 1.0 - np.exp(-k * self.tau)
        else:
            base = 1.0 + k * (self.n - 1.0) * self.tau
            X = 1.0 - (max(base, 1e-10) ** (1.0 / (1.0 - self.n)))
        X = np.clip(X, 0.0, 0.999)
        
        masse_pa6 = feed.imass['PA6_Solid']
        out.imass['PA6_Solid'] -= (masse_pa6 * X)
        out.imass['Caprolactam'] += (masse_pa6 * X * 0.90)
        out.imass['Oligomers'] += (masse_pa6 * X * 0.10)

def simuler_usine_alcoolyse(flux_masses_dictionnaire, type_reacteur, temp_user, pression_user, nom_catalyseur, temps_sejour, nom_solvant, params=None):
    bst.main_flowsheet.clear() 
    bst.settings.set_thermo(chemin_chimique)
    
    ALCOOLS_FACILES = ["2-propanol", "2-butanol", "cyclohexanol", "cyclododecanol", "2-methyl-2-propanol"]
    est_alcool_secondaire = any(alc in str(nom_solvant).lower() for alc in ALCOOLS_FACILES)
    
    soupe_froide = bst.Stream('Melange_PA6_Alcool', **flux_masses_dictionnaire, units='kg/hr', T=25+273.15, P=pression_user * 1000)
    H1_Prechauffage = bst.HXutility('H1_Prechauffage', ins=soupe_froide, outs='Soupe_Chaude', T=temp_user+273.15)
    
    R101 = ReacteurCinetiqueGenerique(ID='R101', ins=H1_Prechauffage.outs[0], outs='Jus_Reactionnel', T_op=temp_user, tau_h=temps_sejour, params=params)
    V1_Valve = bst.IsenthalpicValve('V1_Detente', ins=R101.outs[0], outs='Jus_Detendu', P=101325)
    H2_Flash = bst.HXutility('H2_Flash_Heater', ins=V1_Valve.outs[0], outs='Jus_A_Flasher', T=140+273.15)
    
    Flash_Desalcoolisation = bst.Splitter('Flash_Desalcoolisation', ins=H2_Flash.outs[0], outs=('Vapeur_Alcool', 'Jus_Lourd'), 
                             split={'Isopropanol': 0.999, 'Water': 0.99, 'Caprolactam': 0.01})
    
    C1_Condenseur = bst.HXutility('C1_Condenseur_Alcool', ins=Flash_Desalcoolisation.outs[0], outs='Alcool_Liquide', T=40+273.15)
    Purge_Alcool = bst.Splitter('Purge_Alcool', ins=C1_Condenseur.outs[0], outs=('Alcool_Recycle', 'Purge_Dechets_Legers'), split=0.98)
    C2_Refroidissement = bst.HXutility('C2_Cooler', ins=Flash_Desalcoolisation.outs[1], outs='Jus_Froid', T=80+273.15)
    
    if est_alcool_secondaire:
        Filtre_Solides = bst.Splitter('Filtre_Solides', ins=C2_Refroidissement.outs[0], outs=('Dechets_Solides', 'CPL_Pur'), 
                              split={'PA6_Solid': 1.0, 'Catalyseur_Solide': 1.0, 'Fibre_Verre': 1.0, 'Oligomers': 1.0, 'Caprolactam': 0.01})
        masse_cpl_pur = Filtre_Solides.outs[1].imass['Caprolactam']
        cout_reste_usine = 450000 * (soupe_froide.F_mass / 10000)**0.6 
    else:
        Filtre_Solides = bst.Splitter('Filtre_Solides', ins=C2_Refroidissement.outs[0], outs=('Dechets_Solides', 'Jus_Filtre'), 
                              split={'PA6_Solid': 1.0, 'Catalyseur_Solide': 1.0, 'Fibre_Verre': 1.0, 'Oligomers': 0.0})
        Purification_Finale = bst.Splitter('Purification_Finale', ins=Filtre_Solides.outs[1], outs=('CPL_Impur', 'Fractions_Lourdes'), 
                              split={'Caprolactam': 0.60, 'Oligomers': 0.90})
        masse_cpl_pur = Purification_Finale.outs[0].imass['Caprolactam']
        cout_reste_usine = 2500000 * (soupe_froide.F_mass / 10000)**0.6 

    systeme_alcoolyse = bst.main_flowsheet.create_system('Usine_Alcoolyse')
    systeme_alcoolyse.simulate() 

    try:
        energie_chaude_kJ = soupe_froide.F_mass * 3.0 * max(0.0, temp_user - 25.0) * 1.2 
        volume_reacteur_m3 = (soupe_froide.F_mass / 800.0) * temps_sejour
        cout_equipement_hr = (350000 * (max(volume_reacteur_m3, 0.5) / 5.0)**0.6 + cout_reste_usine) / (10 * 8000)
        masse_appoint_neuf = soupe_froide.imass['Isopropanol'] * 0.02 
        
        gwp_final = ((energie_chaude_kJ/1000)*0.0732 + (energie_chaude_kJ/1000)*0.117 + masse_appoint_neuf*2.50) / masse_cpl_pur if masse_cpl_pur > 0 else 0.0

        return {
            "Cout_Vapeur_hr": (energie_chaude_kJ / 1000.0) * 0.005, "Cout_Eau_Froide_hr": (energie_chaude_kJ / 1000.0) * 0.0005,
            "Eaux_Usees_kg": Purge_Alcool.outs[1].F_mass, "Dechets_Toxiques_kg": Filtre_Solides.outs[0].F_mass,  
            "Cout_Equipement_hr": cout_equipement_hr, "Cout_Solvant_hr": masse_appoint_neuf * (1.60 if est_alcool_secondaire else 0.80),
            "Cout_Reactifs_hr": 0.0, "Cout_Catalyseur_hr": flux_masses_dictionnaire.get('Catalyseur_Solide', 0) * 1.50,
            "Volume_Reacteur_m3": volume_reacteur_m3,
            "Purete_CPL": 95.0, "Type_Usine": "Alcoolyse_Cinetique", "GWP_kgCO2_per_kgCPL": gwp_final, "Production_CPL_kg_hr": masse_cpl_pur
        }
    except Exception as e:
        return {"Type_Usine": "Erreur", "Production_CPL_kg_hr": 0.0}