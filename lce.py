# -----------------------------------------
# Import required Python libraries
# -----------------------------------------
import requests                   # for making HTTP requests to fetch data from IRSA API
import pandas as pd               # for handling tabular data (reading, saving, processing)
from astropy.coordinates import SkyCoord  # for converting RA/Dec into usable formats
import astropy.units as u         # for specifying units like degrees, arcseconds, hours
from io import StringIO           # to treat text response from web as a file (for pandas)

# -----------------------------------------
# Function to query ZTF (Zwicky Transient Facility) lightcurves
# -----------------------------------------
def get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=1.2, band="r", fmt="csv"):
    """
    Query the ZTF IRSA API for lightcurves around a given RA, Dec position.
    Parameters:
        ra_deg (float)  → Right Ascension in degrees
        dec_deg (float) → Declination in degrees
        radius_arcsec   → Search radius around target (arcseconds)
        band (str)      → Which photometric band to fetch ('g', 'r', or 'i')
        fmt (str)       → Data format to request (here CSV)
    """
    base_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"  # API endpoint
    
    # Convert radius from arcseconds to degrees (since API expects degrees)
    radius_deg = radius_arcsec / 3600.0  
    
    # Parameters to send to API
    params = {
        "POS": f"CIRCLE {ra_deg} {dec_deg} {radius_deg}",  # circular search region
        "BANDNAME": band,                                  # photometric filter
        "FORMAT": fmt.lower()                              # response format (csv)
    }
    
    # Make the request to IRSA API
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # raise error if request fails
    
    # Read response text as CSV into a pandas dataframe
    return pd.read_csv(StringIO(response.text), comment="#")

# -----------------------------------------
# Step 1: Load your catalogue of objects
# -----------------------------------------
catalogue = pd.read_csv("magnetic_cataclysmic_variables.txt")  
# (Assumption: file contains columns like "Name", "RA(J2000)", "DEC(J2000)")

# List to hold lightcurves of all objects
all_lightcurves = []

# -----------------------------------------
# Step 2: Loop through every object in catalogue
# -----------------------------------------
for _, row in catalogue.iterrows():   # go row by row in the catalogue
    name = row["Name"]                # object name
    ra_str, dec_str = row["RA(J2000)"], row["DEC(J2000)"]  # RA and Dec as strings
    
    try:
        # Convert RA/Dec strings into decimal degrees using astropy
        coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg
        
        # Print status message
        print(f"Querying ZTF for {name} at RA={ra_deg:.4f}, Dec={dec_deg:.4f}")
        
        # Fetch r-band lightcurve for this object
        df = get_ztf_lightcurve(ra_deg, dec_deg, radius_arcsec=2, band="r")
        
        # If data is returned, attach object name and store
        if not df.empty:
            df["object"] = name
            all_lightcurves.append(df)
            print(f"  ✅ Got {len(df)} points for {name}")
        else:
            print(f"  ⚠ No data returned for {name}")
    
    except Exception as e:
        # Handle any failure (e.g. bad coordinates, API failure, etc.)
        print(f"  ❌ Failed for {name}: {e}")

# -----------------------------------------
# Step 3: Save results into a single CSV file
# -----------------------------------------
if all_lightcurves:
    # Merge all individual lightcurves into one big dataframe
    merged = pd.concat(all_lightcurves, ignore_index=True)
    
    # Save to file
    merged.to_csv("lightcurves_all_r.csv", index=False)
    
    print(f"\n✅ Saved r-band lightcurves for {len(catalogue)} stars ({len(merged)} total points)")
else:
    print("\n❌ No lightcurves fetched.")

