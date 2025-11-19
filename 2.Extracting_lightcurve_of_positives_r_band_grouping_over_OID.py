# ==========================================================
# Import required Python libraries
# ==========================================================
import requests                         # for HTTP requests to fetch data from IRSA API
import pandas as pd                     # for handling tabular data
from astropy.coordinates import SkyCoord  # for converting RA/Dec into degrees
import astropy.units as u               # for specifying units
from io import StringIO                 # to treat text response as a file (for pandas)

# ==========================================================
# Function to query ZTF (Zwicky Transient Facility) lightcurves
# ==========================================================
def get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=1.2, band="r", fmt="csv"):
    """
    Query the ZTF IRSA API for lightcurves around a given RA, Dec position.

    Parameters:
        ra_deg (float)  → Right Ascension in degrees
        dec_deg (float) → Declination in degrees
        radius_arcsec   → Search radius around target (arcseconds)
        band (str)      → Photometric band to fetch ('g', 'r', 'i')
        fmt (str)       → Response format (default: CSV)
    """
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
    
    # Convert radius from arcsec → degrees
    radius_deg = radius_arcsec / 3600.0
    
    # API query parameters
    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",
        "BANDNAME": band,
        "FORMAT": fmt.lower()
    }
    
    # Fetch data
    response = requests.get(base_url, params=params)
    response.raise_for_status()
    
    # Return as pandas DataFrame
    return pd.read_csv(StringIO(response.text), comment="#")

# ==========================================================
# Step 1: Load catalogue of objects
# ==========================================================
catalogue = pd.read_csv("magnetic_cataclysmic_variables.txt")
# Assumption: catalogue has columns → "Name", "RA(J2000)", "DEC(J2000)"

all_lightcurves = []   # container for results

# ==========================================================
# Step 2: Loop through catalogue and fetch lightcurves
# ==========================================================
for _, row in catalogue.iterrows():   # limit to first 5 for testing
    name = row["Name"]
    ra_str, dec_str = row["RA(J2000)"], row["DEC(J2000)"]
    
    try:
        # Convert RA/Dec to decimal degrees
        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg
        
        print(f"Querying ZTF for {name} at RA={ra_deg:.4f}, Dec={dec_deg:.4f}")
        
        # Get r-band lightcurve
        df = get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=1.2, band="r")
        
        if not df.empty:
            df["object"] = name
            all_lightcurves.append(df)
            print(f"  ✅ Got {len(df)} points for {name}")
        else:
            print(f"  ⚠ No data returned for {name}")
    
    except Exception as e:
        print(f"  ❌ Failed for {name}: {e}")

# ==========================================================
# Step 3: Save raw merged results
# ==========================================================
if all_lightcurves:
    merged = pd.concat(all_lightcurves, ignore_index=True)
    merged.to_csv("lightcurves_fetched_from_CV_catalogue.csv", index=False)
    
    print(f"\n✅ Saved r-band lightcurves for {len(catalogue)} stars "
          f"({len(merged)} total points)")
    
    # ==========================================================
    # Step 4: Group by object into compact form (CSV-friendly)
    # ==========================================================
    data = pd.read_csv("lightcurves_fetched_from_CV_catalogue.csv")
    grouped = data.groupby("object", sort=False, as_index=False).agg(list)
    
    # Convert lists → single-line strings (joined by ;)
    for col in grouped.columns:
        if col != "object":
            grouped[col] = grouped[col].apply(lambda x: ";".join(map(str, x)))
    
    grouped.to_csv("lightcurves_grouped_r_clean.csv", index=False)
    print("✅ Grouped output saved to 'lightcurves_grouped_r_clean.csv'")
    
else:
    print("\n❌ No lightcurves fetched.")
