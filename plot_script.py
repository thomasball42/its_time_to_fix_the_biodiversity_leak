# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:49:17 2024

@author: Thomas Ball
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

import _funcs_data

onedrive_path = "C:\\Users\\Thomas Ball\\OneDrive - University of Cambridge"
# onedrive_path = "E:\\OneDrive\\OneDrive - University of Cambridge"

rpath = "D:\\Food_v0\\all_results"  # From NatFood (in revision) 10.33774/coe-2024-fl5fk

datpath = "data\\country_bd_items_weights.csv" # From NatFood (in revision) 10.33774/coe-2024-fl5fk

area_km2 = 1000 #km2 

mix = ["Barley", "Rapeseed", "Wheat"]
mix_yields = [66893, 37390, 85904] # 100g /ha from FAOSTAT
mix_yields_kg_km2 = np.array(mix_yields) * 10 # kg/km2

pdat = pd.read_csv(os.path.join("data", "Production_Crops_Livestock_E_All_Data_(Normalized).csv"),
                   encoding = "latin-1") # FAOSTAT all production data

pdat = pdat[pdat.Year==pdat.Year.unique().max()]
ydat = pdat[pdat.Element.str.contains("Yield")]

if not os.path.isfile(os.path.join("data", "cached_results.csv")) or True:
    
    df = pd.read_csv(datpath, index_col = 0)
     
    ukcode = 229
    
    ukd = pd.read_csv(os.path.join(rpath, "gbr", "df_domestic.csv"))
    uko = pd.read_csv(os.path.join(rpath, "gbr", "df_offshore.csv"))
    
    uko.columns = ["Item"] + uko.columns[1:].to_list()
    
    dfmix = uko[uko["Item"].str.contains("|".join(mix), case = False)]
    
    ukxdf = pd.read_csv(os.path.join(rpath, "gbr", "xdf.csv"), index_col = 0)
    
    dfmix.loc[:, "yields_kg_km2"] = mix_yields_kg_km2
    
    dfmix.loc[:, "disp_prod_kg"] = (1/len(mix)) * area_km2 * dfmix.yields_kg_km2
    
    ldf = pd.DataFrame()
    
    for crop in dfmix.Item:
        
        cdat = ukxdf[(ukxdf.Item == crop)&(ukxdf.Consumer_Country_Code==ukcode)&(ukxdf.Producer_Country_Code!=ukcode)]
        cdat.loc[:, "import_weight"] = cdat.provenance / cdat.provenance.sum() 
        cdat.loc[:, "idisp_prod"] = cdat.loc[:, "import_weight"] * dfmix[dfmix.Item==crop].disp_prod_kg.squeeze()
    
        cropbd = df[df.Item==crop]
        
        for c in cdat.Country_ISO:
            
            bdintensity = cropbd[cropbd.a==c].bd.squeeze()
            if type(bdintensity)==pd.core.series.Series:
                if len(bdintensity) ==0:
                    bdintensity = np.sum(cropbd.bd * cropbd.w) / cropbd.w.sum()
                
            cdat.loc[cdat.Country_ISO==c, "bd_leakage"] = cdat.loc[cdat.Country_ISO==c].idisp_prod.squeeze() * bdintensity
    
        ldf = pd.concat([ldf, cdat])
    
    uk_benefit = df[(df.Item.str.contains("|".join(mix), case=False))&(df.a=="GBR")].bd.mean() * dfmix.disp_prod_kg.sum()
    uk_market_leakage = ldf.bd_leakage.sum()
    
    #%%
    ciso3 = "BRA"
    
    comm_of_interest = "Soy"
    
    ccodes = pd.read_csv(
        os.path.join("data", "country_codes.csv"),
        encoding = "latin-1")
    coi_code = ccodes[ccodes.ISO3==ciso3].FAOSTAT.squeeze()
    coi_exports = pd.DataFrame()
    comm_bd = df[df.Item.str.contains(comm_of_interest,case=False)]
    
    comm_yield_hgha = ydat[(ydat["Area Code"]==coi_code)&(ydat.Item.str.contains(comm_of_interest, case=False))].Value.squeeze()
    comm_yield_kg_km2 = comm_yield_hgha * 10
    
    cdisplaced_production = area_km2 * comm_yield_kg_km2
    
    #%%
    for p in os.walk(rpath):
        path = p[0]
        if path==rpath:
            continue 
        try:
            xdf = pd.read_csv(os.path.join(path, "xdf.csv"), index_col=0)
            cof = xdf[xdf.Item.str.contains(comm_of_interest, case = False)]
            ccof = cof[cof.Producer_Country_Code==coi_code]
            current_country_iso = os.path.split(path)[-1].upper()
            if len(ccof) > 0:
                import_ratio = cof.groupby("Country_ISO").provenance.sum() / cof.provenance.sum()
                cof.loc[:, "export_val"] = cof[cof.Producer_Country_Code==coi_code].provenance.sum()
                for c in cof.Country_ISO.unique():
                    cof.loc[cof.Country_ISO==c, "import_ratio"] = import_ratio[c]       
                cof.loc[:, "Consumer_ISO3"] = current_country_iso
                coi_exports = pd.concat([coi_exports, cof])
        except FileNotFoundError:
            continue
    coi_exports.loc[:, "export_weight"] = coi_exports.export_val / coi_exports.export_val.unique().sum()
    
    # %%
    dpdf = pd.DataFrame()
    for icountry in coi_exports.Consumer_ISO3.unique():    
        cdat = coi_exports[(coi_exports.Consumer_ISO3==icountry)&~(coi_exports.Producer_Country_Code==coi_code)]
        cdat.loc[:, "idisp_prod"] = cdisplaced_production * cdat.export_weight * cdat.import_ratio / cdat.import_ratio.sum()
        dpdf = pd.concat([dpdf, cdat])
        
    prod_sums = dpdf.groupby(["Producer_Country_Code", "Country_ISO"]).idisp_prod.sum().reset_index()
    for c in prod_sums.Country_ISO:
        
        prod_sums.loc[prod_sums.Country_ISO==c, "bd_val"] = comm_bd[comm_bd.a==c].bd.squeeze()
        prod_sums.loc[prod_sums.Country_ISO==c, "bd_leakage"] = prod_sums.loc[prod_sums.Country_ISO==c].idisp_prod.squeeze() * comm_bd[comm_bd.a==c].bd.squeeze()
        
    col_bd = comm_bd[comm_bd.a==ciso3].bd.squeeze()
    col_benefit = comm_bd[comm_bd.a==ciso3].bd.squeeze() * cdisplaced_production
    col_market_leakage = np.nansum(prod_sums.bd_leakage)
    
    #%%
    ad = np.array([uk_benefit, uk_market_leakage] + [col_benefit, col_market_leakage])
    labels = ["Domestic gain", "Market leakage loss"]
    dat = pd.DataFrame([[uk_benefit, col_benefit], [uk_market_leakage, col_market_leakage]], labels, columns = ["UK", "COL"]).to_csv(
        os.path.join(
            "data", "cached_results.csv"))
else:
    dat = pd.read_csv(os.path.join( 
        "data", "cached_results.csv"), index_col = 0)
    
#%% plot
fig, saxs = plt.subplots(2, 3)

a = 0.8
#UK

axs = saxs[:, 0]

bd_km2 = ldf.bd_opp_cost_m2 * 1000000

ax = axs[0]
ax.bar(0, uk_benefit, color = "b", alpha = a)
ax.bar(1, -uk_market_leakage, color = "orange", alpha = a)
xlim = ax.get_xlim()
ax.hlines(0, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)
ax.set_xticks([])
ax.set_ylabel("Avoided extinctions ($\\Delta$E)")
prod_sums = prod_sums 
ldf = ldf # UK
ax = axs[1]
ax.bar(0, col_benefit, color = "b", alpha = a, label = "Domestic gain")
ax.bar(1, -col_market_leakage, color = "orange", alpha = a, label = "Market leakage")
xlim = ax.get_xlim()
ax.hlines(0, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)
ax.set_xticks([])
ax.set_ylabel("Avoided extinctions ($\\Delta$E)")
S = 0.5
axs = saxs[:, 1]
a = 0.8

# UK 
ax = axs[0]

ax.set_title("Mean yield T/km$^2$", size = 10)
yields = ldf.FAO_yield_kgm2
yields_T_km2 = ldf.FAO_yield_kgm2 * 1000000 / 1000
mix_yields_T_km2 = mix_yields_kg_km2 / 1000

ax.boxplot([mix_yields_T_km2.mean()], positions = [S], whis = None, conf_intervals=None, vert = False)
ax.boxplot(yields_T_km2, positions = [0], vert = False, showfliers = False)
ax.set_yticks([])
ax.set_xlim(200, ax.get_xlim()[1])
xlim = ax.get_xlim()
ax.hlines(S/2, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)
ax = axs[1]
ax.set_xlabel("Mean yield T/km$^2$")
yields = dpdf.FAO_yield_kgm2
yields_T_km2 = dpdf.FAO_yield_kgm2 * 1000000 / 1000
coff_T_km2 =  comm_yield_kg_km2 / 1000
ax.boxplot([coff_T_km2], positions = [S], whis = None, conf_intervals=None, vert = False)
ax.boxplot(yields_T_km2, positions = [0], vert = False, showfliers = False)
ax.set_yticks([])
xlim = ax.get_xlim()
ax.hlines(S/2, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)

# BDS
axs = saxs[:, 2]
a = 0.8

# UK
ax = axs[0]
ax.set_yticks([])

UK_bd_km2 = df[(df.Item.str.contains("|".join(mix), case=False))&(df.a=="GBR")].bd * dfmix.yields_kg_km2.mean() / 1000000
ax.boxplot([UK_bd_km2 * 100000], positions = [S], whis = None, conf_intervals=None, vert = False)
ax.set_title("10$^{-5}$$\\Delta$E/km$^2$", size = 10)
ukml_km2 = ldf.groupby("Producer_Country_Code").bd_opp_cost_m2.mean().dropna() * 1000000

ax.boxplot(ukml_km2 * 100000, positions = [0], vert = False, showfliers = False) #whis = (0,100))
ax.set_yticks([])
xlim = ax.get_xlim()
ax.hlines(S/2, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)

#COL
ax = axs[1]
bds = dpdf.groupby("Producer_Country_Code").bd_opp_cost_m2.mean().dropna()
col_bd_km2 = col_bd * comm_yield_kg_km2
bds_km2 = bds * 1000000
ax.set_xlabel("10$^{-5}$$\\Delta$E/km$^2$")

ax.boxplot([col_bd_km2 * 100000], positions = [S], whis = None, conf_intervals=None, vert = False)
ax.boxplot(bds_km2 * 100000, positions = [0], vert = False, whis = (0,100))
xlim = ax.get_xlim()
ax.hlines(S/2, *xlim, color = "k", linewidth = 1)
ax.set_xlim(xlim)
ax.set_yticks([])

fig.tight_layout()