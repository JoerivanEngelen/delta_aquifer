import pandas as pd
import geopandas as gpd
from pathlib import Path
import os


def calc_area_sector(l, phi):
    return phi / 2 * l ** 2


#%% Path management
script_dir = Path(os.path.dirname(os.path.realpath(__file__)))
data_dir = script_dir / ".." / ".." / "data" / "validation"

#%% Process cross-sections
colnames = ["delta", "aqf_nr", "l_start", "l_end", "D", "phi"]
df_cross = pd.DataFrame(columns=colnames)

## Fig 2 in Sarker et al. 2018
name = "Ganges-Brahmaputra"
scale = 100_000 / 4.1
extent = scale * 21  # Total length figure
l_fw = scale * 15.3  # length towards fresh-salt interface aqf 1 & 2
l_hinge = scale * 10  # length towards hinge zone

df_cross.loc[0] = [name, 1, 0, l_hinge, 37.5, 0.31]
df_cross.loc[1] = [name, 1, l_hinge, l_fw, 80, 0.31]
df_cross.loc[2] = [name, 2, l_hinge, l_fw, 75, 0.31]
df_cross.loc[3] = [name, 3, l_hinge, extent, 75, 0.31]

## Compute area and volume based on cross-section and phi
df_cross["A"] = calc_area_sector(df_cross["l_end"], df_cross["phi"]) - calc_area_sector(
    df_cross["l_start"], df_cross["phi"]
)

df_cross["V_fw,on"] = df_cross["A"] * df_cross["D"]

#%% Process map plots
shp_paths = list((data_dir / "shapes").glob("*/*.shp"))
deltas = [i.parent.name for i in shp_paths]

# Aquifer thicknesses estimated from cross-sections in the sources
# See sources.md
thicknesses = [60, 30, 50, 30]

gdf = pd.concat([gpd.read_file(shp_path) for shp_path in shp_paths])
gdf = gdf.reset_index(drop=True)

# Project to equal area projection and compute area of polygon
gdf["area"] = gdf.to_crs({"proj": "cea"}).area

# Add delta
gdf["delta"] = deltas

# Add thicknesses
gdf["D"] = thicknesses

# Compute volume fresh water
gdf["V_fw,on"] = gdf["area"] * gdf["D"]

V_fw_on = gdf.groupby("delta").sum()["V_fw,on"]

# %% Add Ganges
V_fw_on.loc["Ganges-Brahmaputra"] = df_cross["V_fw,on"].sum()
# %% Write to csv
V_fw_on.to_csv(data_dir / "FW_vol_estimations.csv")
