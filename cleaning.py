import pandas as pd
import json
import os

print("--- 1. CHARGEMENT ET NETTOYAGE DES DONNÉES (TRAITEMENT) ---")
dossier_travail = r"c:\Users\mathi\Desktop\Stage_Nylon\python\process"
fichier_source = os.path.join(dossier_travail, 'hydrolyse.csv')

try:
    df = pd.read_csv(fichier_source, sep=',')
    if len(df.columns) < 2:
        df = pd.read_csv(fichier_source, sep=';')
    print(f"✅ Fichier initial 'hydrolyse.csv' chargé avec succès !")
except Exception as e:
    print(f"❌ Erreur de lecture : {e}")
    exit()

df.columns = df.columns.str.strip()

# Correction des virgules décimales
for col in df.columns:
    if df[col].dtype == 'object' and df[col].astype(str).str.contains(',').any():
        df[col] = df[col].astype(str).str.replace(',', '.', regex=False)
    if any(kw in col for kw in ['ºC', 'KPa', 'Time', 'wt%', 'mL/g', 'Yield']):
        df[col] = pd.to_numeric(df[col], errors='coerce')

df['Additive_type'] = df['Additive_type'].fillna('None')
df['Additive_value(wt%)'] = df['Additive_value(wt%)'].fillna(0)
df = df.fillna(0)

print("--- 2. GÉNÉRATION DE L'ARBRE CHIMIQUE (JSON) ---")
dictionnaire_arbre = {}
if 'Process_Type' in df.columns and 'Cat1_type' in df.columns:
    for proc in df['Process_Type'].unique():
        df_sub = df[df['Process_Type'] == proc]
        dictionnaire_arbre[proc] = df_sub['Cat1_type'].unique().tolist()
    
    chemin_json = os.path.join(dossier_travail, "structure_arbre_chimique.json")
    with open(chemin_json, "w", encoding="utf-8") as f:
        json.dump(dictionnaire_arbre, f, indent=4, ensure_ascii=False)

print("--- 3. APPLICATION DU ONE-HOT ENCODING GLOBAL ---")
if 'Source' in df.columns:
    df['Src_Chen'] = (df['Source'].str.strip() == 'Chen').astype(int)
    df = df.drop(columns=['Source'])

# Correction ici pour éviter le warning de Pandas 3/4 en Python 3.14
colonnes_textuelles = df.select_dtypes(include=['object', 'string']).columns.tolist()
df_encoded = pd.get_dummies(df, columns=colonnes_textuelles, drop_first=False)

for col in df_encoded.columns:
    if df_encoded[col].dtype == 'bool':
        df_encoded[col] = df_encoded[col].astype(int)

print("--- 4. EXPORTATION DU FICHIER UNIQUE PAR COMPARTIMENT ---")
chemin_sortie = os.path.join(dossier_travail, 'donnees_preparees.csv')
df_encoded.to_csv(chemin_sortie, sep=';', index=False)

print(f"✅ Terminé ! Nombre de features : {df_encoded.shape[1] - 1}")
print(f"✅ Fichier sauvegardé sous : {chemin_sortie}")