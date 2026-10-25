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
        
        if feed.imass['m_ZrO2'] > 0: 
            cat = feed.imass['m_ZrO2'] / masse_tot * 100
        elif feed.imass['Fe_MOF_Solid'] > 0: 
            cat = feed.imass['Fe_MOF_Solid'] / masse_tot * 100
        
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

def simuler_usine_neutre(flux_masses_dictionnaire, type_reacteur, temp_user, pression_user, nom_catalyseur, temps_sejour, params=None):
    bst.main_flowsheet.clear()
    bst.settings.set_thermo(chemin_chimique)
    
    soupe_froide = bst.Stream('Soupe_Froide', **flux_masses_dictionnaire, units='kg/hr', T=25+273.15, P=pression_user * 1000)
    soupe_prechauffee = soupe_froide.copy()
    soupe_prechauffee.T = 25 + (temp_user - 25) * 0.85 + 273.15
    H1_Trim = bst.HXutility('H1_Trim_Heater', ins=soupe_prechauffee, outs='Soupe_Chaude', T=temp_user+273.15)
    
    R101 = ReacteurCinetiqueGenerique(ID='R101', ins=H1_Trim.outs[0], outs='Jus_Reacteur', T_op=temp_user, tau_h=temps_sejour, params=params)
    C1_Quench = bst.HXutility('C1_Quench_Cooler', ins=R101.outs[0], outs='Soupe_Quench', T=80+273.15)

    Filtre = bst.Splitter('Filtre', ins=C1_Quench.outs[0], outs=('Cat_et_Inertes', 'Jus'), 
                          split={'PA6_Solid': 1.0, 'm_ZrO2': 1.0, 'Fe_MOF_Solid': 1.0, 'Fibre_Verre': 1.0,
                                 'Water': 0.0, 'Caprolactam': 0.0, 'Oligomers': 0.0})
    
    solvant = bst.Stream('Solv', Cyclohexanol=soupe_froide.F_mass * 0.5, units='kg/hr', T=55+273.15)
    Mix = bst.Mixer('Mix', ins=(Filtre.outs[1], solvant), outs=('Melange',))
    
    ELL = bst.Splitter('ELL', ins=Mix.outs[0], outs=('Org', 'Aq'), 
                       split={'Cyclohexanol': 0.99, 'Caprolactam': 0.95, 'Water': 0.05, 'Oligomers': 0.01})
    
    H2_Reboiler = bst.HXutility('H2_Reboiler', ins=ELL.outs[0], outs='Org_Chaud', T=165+273.15)
    Flash = bst.Splitter('Flash', ins=H2_Reboiler.outs[0], outs=('Vapeur_Solvant', 'CPL_Brut'), 
                         split={'Cyclohexanol': 0.99, 'Caprolactam': 0.01, 'Water': 0.99, 'Oligomers': 0.0})
    
    C2_Condenser = bst.HXutility('C2_Condenser', ins=Flash.outs[0], outs='Solvant_Froid', T=55+273.15)
    Purge = bst.Splitter('Purge', ins=C2_Condenser.outs[0], outs=('Solv_Recycle', 'Purge_Dechets'), split=0.98)
    C3_Chiller = bst.HXutility('C3_Chiller', ins=Flash.outs[1], outs='CPL_Froid', T=60+273.15)
    Cristaux = bst.Splitter('Cristaux', ins=C3_Chiller.outs[0], outs=('CPL_Pur', 'Dechets_Liq'), 
                            split={'Caprolactam': 0.98, 'Cyclohexanol': 0.001, 'Water': 0.001})
    
    mon_systeme_complet = bst.main_flowsheet.create_system('Usine_Neutre')
    mon_systeme_complet.simulate()

    try:
        energie_flash = H2_Reboiler.outs[0].H - H2_Reboiler.ins[0].H
        masse_appoint_neuf = max(0.0, (soupe_froide.F_mass * 0.5) - Purge.outs[0].imass['Cyclohexanol'])
        
        masse_entree_totale = soupe_froide.F_mass 
        energie_chaude_totale_kJ = (masse_entree_totale * 4.18 * max(0.0, temp_user - 25.0) * 1.2) + abs(energie_flash)
        
        cout_vapeur_hr = (energie_chaude_totale_kJ / 1000.0) * 0.005
        cout_eau_froide_hr = (energie_chaude_totale_kJ / 1000.0) * 0.0005
        volume_reacteur_m3 = (masse_entree_totale / 950.0) * temps_sejour
        cout_equipement_hr = (350000 * (max(volume_reacteur_m3, 0.5) / 5.0)**0.6 + 2150000 * (masse_entree_totale / 10000)**0.6) / (10 * 8000)
        
        masse_dechets = ELL.outs[1].F_mass + Cristaux.outs[1].F_mass + Purge.outs[1].F_mass + Filtre.outs[0].imass['Fibre_Verre']
        masse_cata_totale = flux_masses_dictionnaire.get('m_ZrO2', 0) + flux_masses_dictionnaire.get('Fe_MOF_Solid', 0)
        cout_cata_hr = (masse_cata_totale * 0.02) * {"Fe- MOF-1": 150.0, "Fe- MOF-2": 150.0, "m-ZrO2": 80.0}.get(nom_catalyseur, 0.0)
        
        masse_cpl_pur = Cristaux.outs[0].imass['Caprolactam']
        purete_finale = (masse_cpl_pur / Cristaux.outs[0].F_mass) * 100.0 if Cristaux.outs[0].F_mass > 0 else 0.0
        
        emissions = (energie_chaude_totale_kJ/1000)*0.0732 + (energie_chaude_totale_kJ/1000)*0.117 + masse_appoint_neuf*2.50 + (masse_cata_totale*0.02)*15.0
        gwp_final = emissions / masse_cpl_pur if masse_cpl_pur > 0 else 0.0

        return {
            "Cout_Vapeur_hr": cout_vapeur_hr, "Cout_Eau_Froide_hr": cout_eau_froide_hr, "Dechets_Totaux_kg": masse_dechets,          
            "Type_Usine": "Hydrolyse Neutre", "GWP_kgCO2_per_kgCPL": gwp_final,"Volume_Reacteur_m3": volume_reacteur_m3,
            "Cout_Equipement_hr": cout_equipement_hr, "Cout_Solvant_hr": masse_appoint_neuf * 1.50, "Production_CPL_kg_hr": masse_cpl_pur,
            "Cout_Reactifs_hr": 0.0, "Cout_Catalyseur_hr": cout_cata_hr, "Purete_CPL": purete_finale,
        }
    except Exception as e:
        return {"Type_Usine": "Erreur", "Production_CPL_kg_hr": 0.0}