import geopandas as gpd

# Load the geojson file
gdf = gpd.read_file("toronto.geojson")

# Reproject to UTM zone 17N (an equal-area-ish local projection)
gdf = gdf.to_crs(26917)

# Compute area in km^2
gdf["area_km2"] = gdf.geometry.area / 1_000_000

# Show a preview
print(gdf[["name", "area_km2"]])

# Save to CSV
gdf[["name", "area_km2"]].to_csv("toronto_areas.csv", index=False)
print("\nSaved to toronto_areas.csv")

