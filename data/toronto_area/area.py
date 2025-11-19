import geopandas as gpd
import pandas as pd

# ---- config ----
INPUT_GEOJSON = "toronto_full_neighborhoods.geojson"
OUTPUT_CSV    = "fixed_neighbourhood_areas.csv"   # rename to .sav if you just want a generic save file

# ---- load geojson ----
gdf = gpd.read_file(INPUT_GEOJSON)

# If CRS is missing, assume WGS84 (long/lat)
if gdf.crs is None:
    gdf = gdf.set_crs(epsg=4326)

# Reproject to an equal-area CRS (Canada Albers) and compute area in km^2
gdf_proj = gdf.to_crs(epsg=3978)  # EPSG:3978 = NAD83 / Canada Albers Equal Area
gdf_proj["area_km2"] = gdf_proj.geometry.area / 1_000_000.0

# ---- pick the neighbourhood number and name fields ----
id_candidates = [
    "CLASSIFICATION_CODE",
    "AREA_SHORT_CODE",
    "AREA_LONG_CODE",
    "AREA_ID",
    "NEIGH_ID",
    "id",
]

name_candidates = [
    "AREA_NAME",
    "NEIGHBOURHOOD",
    "Neighbourhood",
    "NAME",
]

# choose ID field
id_field = None
for c in id_candidates:
    if c in gdf_proj.columns:
        id_field = c
        break

# if there is no ID field, just make one 1..N
if id_field is None:
    gdf_proj["neighbourhood_number"] = range(1, len(gdf_proj) + 1)
    id_field = "neighbourhood_number"

# choose name field
name_field = None
for c in name_candidates:
    if c in gdf_proj.columns:
        name_field = c
        break

if name_field is None:
    raise ValueError("Couldn't find a neighbourhood name field in properties.")

# ---- build output table ----
out = gdf_proj[[id_field, name_field, "area_km2"]].copy()
out.columns = ["neighbourhood_number", "neighbourhood_name", "area_km2"]

# ---- save ----
out.to_csv(OUTPUT_CSV, index=False)

print(f"Saved {len(out)} rows to {OUTPUT_CSV}")
print(out.head())
